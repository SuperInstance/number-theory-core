# number-theory-core

A pure-Rust library implementing **core number theory algorithms** including prime
sieving, Miller-Rabin primality testing, integer factorization, divisor functions,
and the Möbius function with inversion.

## Why This Matters

Number theory is one of the oldest and most beautiful branches of mathematics,
and it underpins all of modern cryptography. RSA encryption relies on the
difficulty of factoring large integers. Diffie-Hellman key exchange requires
finding large primes. Elliptic curve cryptography needs efficient primality
testing. The distribution of prime numbers, divisor sums, and the Möbius function
connect deep problems like the Riemann Hypothesis to computational algorithms.

This library provides clean, efficient, and thoroughly tested implementations
of these fundamental algorithms, suitable for both educational exploration and
as building blocks for cryptographic applications.

## Features

- **Sieve of Eratosthenes** — generate all primes up to n in O(n log log n)
- **Prime counting** — count primes up to n
- **Miller-Rabin primality test** — deterministic for n < 3.3 × 10²⁴
- **Trial division** — factorize by exhaustive trial up to √n
- **Pollard's rho** — efficient factorization for medium-sized composites
- **Full factorization** — trial division + Pollard's rho combined
- **Divisor functions** — σₖ(n), number of divisors d(n), sum of divisors σ(n)
- **Möbius function** — μ(n) for individual values and sieve
- **Möbius inversion** — recover g from f = g ∗ 1
- **Euler's totient** — φ(n) via factorization

## Mathematical Background

### Sieve of Eratosthenes

The ancient algorithm for finding all primes up to n. Mark multiples of each
prime starting from p². Time complexity O(n log log n), which is nearly linear.

The prime counting function π(n) ~ n / ln(n) by the Prime Number Theorem.

### Miller-Rabin Primality Test

A probabilistic (or deterministic for bounded inputs) primality test based on
the following property: if n is prime, then for any 0 < a < n, writing n-1 = 2^r · d:

Either a^d ≡ 1 (mod n), or a^(2^j · d) ≡ -1 (mod n) for some 0 ≤ j < r.

If n is composite, a random base a has at least a 3/4 chance of being a witness
to compositeness. Specific witness sets make the test **deterministic** for all
n up to enormous bounds.

### Integer Factorization

**Trial division**: try dividing by 2, 3, 4, ..., up to √n. Simple but O(√n).

**Pollard's rho**: Uses Floyd's cycle-finding algorithm on the iteration
x → x² + c (mod n) to find a non-trivial factor. Expected time O(n^(1/4)).

For the full factorization, we combine both: trial division for small factors,
then Pollard's rho for the remaining composite parts.

### Divisor Functions

The **divisor function** σₖ(n) = Σ_{d|n} d^k counts or sums divisors:

- **d(n) = σ₀(n)**: number of divisors (e.g., d(12) = 6)
- **σ(n) = σ₁(n)**: sum of divisors (e.g., σ(12) = 28)

A number is **perfect** if σ(n) = 2n (e.g., 6, 28, 496, 8128).

For n = p₁^a₁ · p₂^a₂ · ..., σₖ(n) = ∏ (pᵢ^(k(aᵢ+1)) - 1) / (pᵢ^k - 1).

### Möbius Function

The Möbius function μ(n) is defined as:

- μ(1) = 1
- μ(n) = (-1)^k if n is a product of k distinct primes
- μ(n) = 0 if n has a squared prime factor

**Möbius inversion**: If f(n) = Σ_{d|n} g(d), then g(n) = Σ_{d|n} μ(n/d) · f(d).

This is fundamental in analytic number theory and connects to the Riemann
zeta function: 1/ζ(s) = Σ μ(n)/n^s.

## Usage

```toml
[dependencies]
number-theory-core = "0.1.0"
```

### Prime Sieve

```rust
use number_theory_core::{prime_sieve, prime_count};

let primes = prime_sieve(100);
assert_eq!(primes.len(), 25); // π(100) = 25

assert_eq!(prime_count(1000), 168);
```

### Primality Testing

```rust
use number_theory_core::is_prime;

assert!(is_prime(104729));      // 10000th prime
assert!(!is_prime(104730));     // Composite
assert!(!is_prime(561));        // Carmichael number (still detected!)
```

### Factorization

```rust
use number_theory_core::{trial_division, factorize, pollard_rho};

// Trial division for small numbers
let factors = trial_division(360);
// 360 = 2³ × 3² × 5

// Full factorization (trial + Pollard's rho)
let factors = factorize(30030);
// 30030 = 2 × 3 × 5 × 7 × 11 × 13
```

### Divisor Functions

```rust
use number_theory_core::{divisor_count, divisor_sum, sigma};

assert_eq!(divisor_count(12), 6);   // {1,2,3,4,6,12}
assert_eq!(divisor_sum(28), 56);    // 28 is perfect: σ(28) = 56
assert_eq!(sigma(4, 2), 21);       // 1 + 4 + 16 = 21
```

### Möbius Function

```rust
use number_theory_core::{mobius, mobius_sieve};

assert_eq!(mobius(1), 1);   // By convention
assert_eq!(mobius(6), 1);   // 2×3, two distinct primes → (-1)² = 1
assert_eq!(mobius(4), 0);   // 2², squared factor

// Sieve μ for all n ≤ N
let mu = mobius_sieve(100);
assert_eq!(mu[30], -1); // 2×3×5, three primes → (-1)³ = -1
```

### Euler's Totient

```rust
use number_theory_core::euler_totient;

assert_eq!(euler_totient(12), 4);  // {1, 5, 7, 11}
assert_eq!(euler_totient(7), 6);   // φ(p) = p - 1
```

## API Reference

| Function | Description |
|---|---|
| `prime_sieve(n)` | All primes ≤ n (Sieve of Eratosthenes) |
| `prime_count(n)` | Count of primes ≤ n |
| `is_prime(n)` | Miller-Rabin primality test |
| `trial_division(n)` | Factorize via trial division |
| `pollard_rho(n)` | Non-trivial factor via Pollard's rho |
| `factorize(n)` | Full factorization (trial + rho) |
| `sigma(n, k)` | σₖ(n) = Σ d^k for d \| n |
| `divisor_count(n)` | d(n) = σ₀(n) |
| `divisor_sum(n)` | σ(n) = σ₁(n) |
| `mobius(n)` | Möbius function μ(n) |
| `mobius_sieve(limit)` | μ for all n ≤ limit |
| `mobius_inversion_sum(f, n)` | Möbius inversion |
| `euler_totient(n)` | Euler's φ function |

## Testing

```bash
cargo test
```

## License

MIT License. See [LICENSE](LICENSE) for details.
