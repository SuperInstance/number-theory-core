//! # Core Number Theory Algorithms
//!
//! A pure-Rust library implementing fundamental number theory algorithms including
//! prime sieving, primality testing, integer factorization, divisor functions, and
//! the Möbius function.

use std::collections::HashMap;

/// Sieve of Eratosthenes: returns all primes up to n.
///
/// Time complexity: O(n log log n)
/// Space complexity: O(n)
pub fn prime_sieve(n: usize) -> Vec<u64> {
    if n < 2 {
        return vec![];
    }

    let mut is_prime = vec![true; n + 1];
    is_prime[0] = false;
    is_prime[1] = false;

    let mut i = 2;
    while i * i <= n {
        if is_prime[i] {
            let mut j = i * i;
            while j <= n {
                is_prime[j] = false;
                j += i;
            }
        }
        i += 1;
    }

    (2..=n).filter(|&i| is_prime[i]).map(|i| i as u64).collect()
}

/// Count primes up to n using the sieve.
pub fn prime_count(n: usize) -> usize {
    prime_sieve(n).len()
}

/// Miller-Rabin primality test.
///
/// Deterministic for all n < 3,317,044,064,679,887,385,961,981
/// using specific witness sets.
///
/// For larger n, uses probabilistic testing with the given number of rounds.
pub fn is_prime(n: u64) -> bool {
    if n < 2 { return false; }
    if n < 4 { return true; }
    if n.is_multiple_of(2) { return false; }
    if n.is_multiple_of(3) { return false; }

    // Write n - 1 as 2^r * d
    let mut d = n - 1;
    let mut r = 0u32;
    while d.is_multiple_of(2) {
        d /= 2;
        r += 1;
    }

    // Deterministic witnesses for n < 3.3 × 10^24
    let witnesses: &[u64] = if n < 2_047 {
        &[2]
    } else if n < 1_373_653 {
        &[2, 3]
    } else if n < 9_080_191 {
        &[31, 73]
    } else if n < 25_326_001 {
        &[2, 3, 5]
    } else if n < 3_215_031_751 {
        &[2, 3, 5, 7]
    } else if n < 4_759_123_141 {
        &[2, 7, 61]
    } else if n < 1_122_004_669_633 {
        &[2, 13, 23, 1662803]
    } else if n < 2_152_302_898_747 {
        &[2, 3, 5, 7, 11]
    } else if n < 3_474_749_660_383 {
        &[2, 3, 5, 7, 11, 13]
    } else if n < 341_550_071_728_321 {
        &[2, 3, 5, 7, 11, 13, 17]
    } else {
        &[2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    };

    'witness: for &a in witnesses {
        if a >= n { continue; }

        let mut x = mod_pow_u64(a, d, n);

        if x == 1 || x == n - 1 {
            continue 'witness;
        }

        for _ in 0..r - 1 {
            x = mod_pow_u64(x, 2, n);
            if x == n - 1 {
                continue 'witness;
            }
        }

        return false;
    }

    true
}

/// Modular exponentiation for u64.
fn mod_pow_u64(mut base: u64, mut exp: u64, modulus: u64) -> u64 {
    if modulus == 1 { return 0; }
    let mut result = 1u64;
    base %= modulus;
    while exp > 0 {
        if exp & 1 == 1 {
            result = mod_mul_u64(result, base, modulus);
        }
        exp >>= 1;
        if exp > 0 {
            base = mod_mul_u64(base, base, modulus);
        }
    }
    result
}

/// Modular multiplication avoiding overflow for u64.
fn mod_mul_u64(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

/// Trial division factorization.
///
/// Returns prime factors with multiplicity as a sorted Vec of (prime, exponent) pairs.
pub fn trial_division(mut n: u64) -> Vec<(u64, u64)> {
    if n <= 1 { return vec![]; }

    let mut factors = Vec::new();
    let mut d = 2u64;

    while d * d <= n {
        if n.is_multiple_of(d) {
            let mut count = 0u64;
            while n.is_multiple_of(d) {
                n /= d;
                count += 1;
            }
            factors.push((d, count));
        }
        d += 1;
    }

    if n > 1 {
        factors.push((n, 1));
    }

    factors
}

/// Pollard's rho factorization algorithm.
///
/// Returns a non-trivial factor of n, or None if n is prime.
pub fn pollard_rho(n: u64) -> Option<u64> {
    if n <= 1 { return None; }
    if n.is_multiple_of(2) { return Some(2); }
    if is_prime(n) { return None; }

    // Try small factors first
    for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31] {
        if n.is_multiple_of(p) { return Some(p); }
    }

    let n128 = n as u128;

    for c in 1u64.. {
        let f = |x: u64| -> u64 {
            ((x as u128 * x as u128 + c as u128) % n128) as u64
        };

        let mut x = 2u64;
        let mut y = 2u64;
        let mut d = 1u64;

        while d == 1 {
            x = f(x);
            y = f(f(y));
            let diff = x.abs_diff(y);
            d = gcd_u64(diff, n);
        }

        if d != n {
            return Some(d);
        }
    }

    None
}

/// Full factorization using trial division + Pollard's rho.
pub fn factorize(n: u64) -> Vec<(u64, u64)> {
    if n <= 1 { return vec![]; }

    let mut factors: HashMap<u64, u64> = HashMap::new();
    let mut queue = vec![n];

    while let Some(m) = queue.pop() {
        if m <= 1 { continue; }
        if is_prime(m) {
            *factors.entry(m).or_insert(0) += 1;
        } else {
            // Try trial division first for small factors
            let td = trial_division_small(m);
            if let Some((p, _)) = td.first() {
                if *p < 10000 {
                    for (p, e) in td {
                        *factors.entry(p).or_insert(0) += e;
                    }
                    continue;
                }
            }
            // Use Pollard's rho
            if let Some(d) = pollard_rho(m) {
                queue.push(d);
                queue.push(m / d);
            } else {
                // Shouldn't happen if is_prime is correct
                *factors.entry(m).or_insert(0) += 1;
            }
        }
    }

    let mut result: Vec<(u64, u64)> = factors.into_iter().collect();
    result.sort_by_key(|&(p, _)| p);
    result
}

/// Trial division up to a small bound.
fn trial_division_small(mut n: u64) -> Vec<(u64, u64)> {
    let mut factors = Vec::new();
    let mut d = 2u64;
    while d * d <= n && d < 10000 {
        if n.is_multiple_of(d) {
            let mut count = 0u64;
            while n.is_multiple_of(d) {
                n /= d;
                count += 1;
            }
            factors.push((d, count));
        }
        d += 1;
    }
    if n > 1 {
        factors.push((n, 1));
    }
    factors
}

fn gcd_u64(a: u64, b: u64) -> u64 {
    if b == 0 { a } else { gcd_u64(b, a % b) }
}

/// Sum-of-divisors function σ_k(n).
///
/// σ_k(n) = Σ_{d|n} d^k
///
/// For k=0 this gives the number of divisors.
/// For k=1 this gives the sum of divisors.
pub fn sigma(n: u64, k: u64) -> u64 {
    if n == 0 { return 0; }
    let factors = trial_division(n);
    let mut result = 1u64;

    for &(p, e) in &factors {
        // σ_k(p^e) = (p^(k(e+1)) - 1) / (p^k - 1) for k > 0
        // σ_0(p^e) = e + 1
        if k == 0 {
            result *= e + 1;
        } else {
            let mut sum = 1u64;
            let mut pk = 1u64; // p^k
            for _ in 0..k {
                pk = pk.saturating_mul(p);
            }
            let mut term = 1u64;
            for _ in 0..e {
                term = term.saturating_mul(pk);
                sum = sum.saturating_add(term);
            }
            result = result.saturating_mul(sum);
        }
    }

    result
}

/// Number of divisors d(n) = σ₀(n).
pub fn divisor_count(n: u64) -> u64 {
    sigma(n, 0)
}

/// Sum of divisors σ(n) = σ₁(n).
pub fn divisor_sum(n: u64) -> u64 {
    sigma(n, 1)
}

/// Möbius function μ(n).
///
/// Returns:
/// - 0 if n has a squared prime factor
/// - 1 if n is a product of an even number of distinct primes
/// - -1 if n is a product of an odd number of distinct primes
pub fn mobius(n: u64) -> i64 {
    if n == 0 { return 0; }
    if n == 1 { return 1; }

    let factors = trial_division(n);
    for &(_, e) in &factors {
        if e > 1 {
            return 0;
        }
    }

    if factors.len().is_multiple_of(2) { 1 } else { -1 }
}

/// Möbius function values for all n in [1, limit].
///
/// Returns a vector where mu[n] = μ(n).
pub fn mobius_sieve(limit: usize) -> Vec<i64> {
    if limit == 0 { return vec![]; }

    let mut mu = vec![0i64; limit + 1];
    mu[1] = 1;

    let mut is_prime = vec![true; limit + 1];
    let mut primes = Vec::new();

    for i in 2..=limit {
        if is_prime[i] {
            primes.push(i);
            mu[i] = -1;
        }

        for &p in &primes {
            if i * p > limit { break; }
            is_prime[i * p] = false;

            if i % p == 0 {
                mu[i * p] = 0;
                break;
            } else {
                mu[i * p] = -mu[i];
            }
        }
    }

    mu
}

/// Möbius inversion: given f where f(n) = Σ_{d|n} g(d),
/// computes g(n) = Σ_{d|n} μ(n/d) · f(d).
///
/// For simplicity, this implements the summatory version:
/// G(x) = Σ_{n≤x} g(n) = Σ_{n≤x} μ(n) · F(x/n)
/// where F(x) = Σ_{n≤x} f(n).
pub fn mobius_inversion_sum(f: &[i64], n: usize) -> i64 {
    let mu = mobius_sieve(n);
    let mut result = 0i64;

    for d in 1..=n {
        result += mu[d] * f.get(n / d).copied().unwrap_or(0);
    }

    result
}

/// Euler's totient function φ(n) via factorization.
pub fn euler_totient(n: u64) -> u64 {
    if n <= 1 { return n; }
    let factors = trial_division(n);
    let mut result = n;
    for &(p, _) in &factors {
        result = result / p * (p - 1);
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_prime_sieve_small() {
        let primes = prime_sieve(30);
        assert_eq!(primes, vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29]);
    }

    #[test]
    fn test_prime_sieve_empty() {
        assert_eq!(prime_sieve(1), vec![]);
        assert_eq!(prime_sieve(0), vec![]);
    }

    #[test]
    fn test_prime_count() {
        assert_eq!(prime_count(100), 25);
        assert_eq!(prime_count(1000), 168);
    }

    #[test]
    fn test_is_prime_small() {
        assert!(!is_prime(0));
        assert!(!is_prime(1));
        assert!(is_prime(2));
        assert!(is_prime(3));
        assert!(!is_prime(4));
        assert!(is_prime(5));
        assert!(is_prime(7));
        assert!(!is_prime(9));
    }

    #[test]
    fn test_is_prime_medium() {
        assert!(is_prime(104729)); // 10000th prime
        assert!(!is_prime(104730));
        assert!(is_prime(999999937)); // Large prime
    }

    #[test]
    fn test_is_prime_carmichael() {
        // Carmichael numbers are composite but Fermat pseudoprimes
        assert!(!is_prime(561));   // 3 × 11 × 17
        assert!(!is_prime(1105));  // 5 × 13 × 17
        assert!(!is_prime(1729));  // 7 × 13 × 19
    }

    #[test]
    fn test_trial_division() {
        let factors = trial_division(12);
        assert_eq!(factors, vec![(2, 2), (3, 1)]); // 12 = 2² × 3
    }

    #[test]
    fn test_trial_division_prime() {
        let factors = trial_division(13);
        assert_eq!(factors, vec![(13, 1)]);
    }

    #[test]
    fn test_trial_division_one() {
        assert_eq!(trial_division(1), vec![]);
    }

    #[test]
    fn test_pollard_rho_composite() {
        let n = 15u64; // 3 × 5
        let factor = pollard_rho(n).unwrap();
        assert!(n % factor == 0 && factor > 1 && factor < n);
    }

    #[test]
    fn test_pollard_rho_prime() {
        assert!(pollard_rho(13).is_none());
    }

    #[test]
    fn test_factorize() {
        let factors = factorize(360);
        // 360 = 2³ × 3² × 5
        assert_eq!(factors, vec![(2, 3), (3, 2), (5, 1)]);
    }

    #[test]
    fn test_divisor_count() {
        // d(12) = 6: divisors are {1, 2, 3, 4, 6, 12}
        assert_eq!(divisor_count(12), 6);
        assert_eq!(divisor_count(1), 1);
        assert_eq!(divisor_count(28), 6); // Perfect number
    }

    #[test]
    fn test_divisor_sum() {
        // σ(12) = 1+2+3+4+6+12 = 28
        assert_eq!(divisor_sum(12), 28);
        // σ(6) = 1+2+3+6 = 12 (perfect number)
        assert_eq!(divisor_sum(6), 12);
    }

    #[test]
    fn test_perfect_number() {
        // 28 is perfect: σ(28) = 56 = 2×28, so σ(28) - 28 = 28
        assert_eq!(divisor_sum(28), 56);
    }

    #[test]
    fn test_sigma_k() {
        // σ₂(4) = 1 + 4 + 16 = 21
        assert_eq!(sigma(4, 2), 21);
    }

    #[test]
    fn test_mobius() {
        assert_eq!(mobius(1), 1);   // μ(1) = 1
        assert_eq!(mobius(2), -1);  // Single prime
        assert_eq!(mobius(4), 0);   // 2², squared factor
        assert_eq!(mobius(6), 1);   // 2×3, two primes
        assert_eq!(mobius(30), -1); // 2×3×5, three primes
        assert_eq!(mobius(12), 0);  // 2²×3, squared factor
    }

    #[test]
    fn test_mobius_sieve() {
        let mu = mobius_sieve(10);
        assert_eq!(mu[1], 1);
        assert_eq!(mu[2], -1);
        assert_eq!(mu[3], -1);
        assert_eq!(mu[4], 0);
        assert_eq!(mu[5], -1);
        assert_eq!(mu[6], 1);
        assert_eq!(mu[7], -1);
        assert_eq!(mu[8], 0);
        assert_eq!(mu[9], 0);
        assert_eq!(mu[10], 1);
    }

    #[test]
    fn test_euler_totient() {
        assert_eq!(euler_totient(1), 1);
        assert_eq!(euler_totient(12), 4); // {1, 5, 7, 11}
        assert_eq!(euler_totient(7), 6);  // Prime
    }

    #[test]
    fn test_factorize_large_composite() {
        let n = 2 * 3 * 5 * 7 * 11 * 13; // 30030
        let factors = factorize(n);
        assert_eq!(factors, vec![(2, 1), (3, 1), (5, 1), (7, 1), (11, 1), (13, 1)]);
    }
}
