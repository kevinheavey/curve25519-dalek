//! Arithmetic mod \\(2\^{252} + 27742317777372353535851937790883648493\\)
//! with five \\(52\\)-bit unsigned limbs.
//!
//! \\(51\\)-bit limbs would cover the desired bit range (\\(253\\)
//! bits), but isn't large enough to reduce a \\(512\\)-bit number with
//! Montgomery multiplication, so \\(52\\) bits is used instead.  To see
//! that this is safe for intermediate results, note that the largest
//! limb in a \\(5\times 5\\) product of \\(52\\)-bit limbs will be
//!
//! ```text
//! (0xfffffffffffff^2) * 5 = 0x4ffffffffffff60000000000005 (107 bits).
//! ```

use core::ops::{Index, IndexMut};
use crate::constants;

/// The `Scalar52` struct represents an element in
/// \\(\mathbb Z / \ell \mathbb Z\\) as 5 \\(52\\)-bit limbs.
#[derive(Copy, Clone)]
pub(crate) struct Scalar52(pub [u64; 5]);

impl Index<usize> for Scalar52 {
    type Output = u64;
    fn index(&self, _index: usize) -> &u64 {
        &(self.0[_index])
    }
}

impl IndexMut<usize> for Scalar52 {
    fn index_mut(&mut self, _index: usize) -> &mut u64 {
        &mut (self.0[_index])
    }
}

/// u64 * u64 = u128 multiply helper
#[inline(always)]
fn m(x: u64, y: u64) -> u128 {
    (x as u128) * (y as u128)
}

impl Scalar52 {
    /// The scalar \\( 0 \\).
    pub(crate) const ZERO: Scalar52 = Scalar52([0, 0, 0, 0, 0]);

    /// Compute `a - b` (mod l)
    pub(crate) fn sub(a: &Scalar52, b: &Scalar52) -> Scalar52 {
        let mut difference = Scalar52::ZERO;
        let mask = (1u64 << 52) - 1;

        // a - b
        let mut borrow: u64 = 0;
        for i in 0..5 {
            borrow = a[i].wrapping_sub(b[i] + (borrow >> 63));
            difference[i] = borrow & mask;
        }

        // conditionally add l if the difference is negative
        let underflow_mask = ((borrow >> 63) ^ 1).wrapping_sub(1);
        let mut carry: u64 = 0;
        for i in 0..5 {
            carry = (carry >> 52) + difference[i] + (constants::L[i] & underflow_mask);
            difference[i] = carry & mask;
        }

        difference
    }

    /// Compute `a * b`
    #[inline(always)]
    #[rustfmt::skip] // keep alignment of z[*] calculations
    pub (crate) fn mul_internal(a: &Scalar52, b: &Scalar52) -> [u128; 9] {
        let mut z = [0u128; 9];

        z[0] = m(a[0], b[0]);
        z[1] = m(a[0], b[1]) + m(a[1], b[0]);
        z[2] = m(a[0], b[2]) + m(a[1], b[1]) + m(a[2], b[0]);
        z[3] = m(a[0], b[3]) + m(a[1], b[2]) + m(a[2], b[1]) + m(a[3], b[0]);
        z[4] = m(a[0], b[4]) + m(a[1], b[3]) + m(a[2], b[2]) + m(a[3], b[1]) + m(a[4], b[0]);
        z[5] =                 m(a[1], b[4]) + m(a[2], b[3]) + m(a[3], b[2]) + m(a[4], b[1]);
        z[6] =                                 m(a[2], b[4]) + m(a[3], b[3]) + m(a[4], b[2]);
        z[7] =                                                 m(a[3], b[4]) + m(a[4], b[3]);
        z[8] =                                                                 m(a[4], b[4]);

        z
    }

    /// Compute `a^2`
    #[inline(always)]
    #[rustfmt::skip] // keep alignment of return calculations
    fn square_internal(a: &Scalar52) -> [u128; 9] {
        let aa = [
            a[0] * 2,
            a[1] * 2,
            a[2] * 2,
            a[3] * 2,
        ];

        [
            m( a[0], a[0]),
            m(aa[0], a[1]),
            m(aa[0], a[2]) + m( a[1], a[1]),
            m(aa[0], a[3]) + m(aa[1], a[2]),
            m(aa[0], a[4]) + m(aa[1], a[3]) + m( a[2], a[2]),
                             m(aa[1], a[4]) + m(aa[2], a[3]),
                                              m(aa[2], a[4]) + m( a[3], a[3]),
                                                               m(aa[3], a[4]),
                                                                                m(a[4], a[4])
        ]
    }

    /// Compute `limbs/R` (mod l), where R is the Montgomery modulus 2^260
    #[inline(always)]
    #[rustfmt::skip] // keep alignment of n* and r* calculations
    pub (crate) fn montgomery_reduce(limbs: &[u128; 9]) -> Scalar52 {

        #[inline(always)]
        fn part1(sum: u128) -> (u128, u64) {
            let p = (sum as u64).wrapping_mul(constants::LFACTOR) & ((1u64 << 52) - 1);
            ((sum + m(p, constants::L[0])) >> 52, p)
        }

        #[inline(always)]
        fn part2(sum: u128) -> (u128, u64) {
            let w = (sum as u64) & ((1u64 << 52) - 1);
            (sum >> 52, w)
        }

        // note: l[3] is zero, so its multiples can be skipped
        let l = &constants::L;

        // the first half computes the Montgomery adjustment factor n, and begins adding n*l to make limbs divisible by R
        let (carry, n0) = part1(        limbs[0]);
        let (carry, n1) = part1(carry + limbs[1] + m(n0, l[1]));
        let (carry, n2) = part1(carry + limbs[2] + m(n0, l[2]) + m(n1, l[1]));
        let (carry, n3) = part1(carry + limbs[3]               + m(n1, l[2]) + m(n2, l[1]));
        let (carry, n4) = part1(carry + limbs[4] + m(n0, l[4])               + m(n2, l[2]) + m(n3, l[1]));

        // limbs is divisible by R now, so we can divide by R by simply storing the upper half as the result
        let (carry, r0) = part2(carry + limbs[5]               + m(n1, l[4])               + m(n3, l[2])   + m(n4, l[1]));
        let (carry, r1) = part2(carry + limbs[6]                             + m(n2,l[4])                  + m(n4, l[2]));
        let (carry, r2) = part2(carry + limbs[7]                                           + m(n3, l[4])                );
        let (carry, r3) = part2(carry + limbs[8]                                                           + m(n4, l[4]));
        let         r4 = carry as u64;

        // result may be >= l, so attempt to subtract l
        Scalar52::sub(&Scalar52([r0, r1, r2, r3, r4]), l)
    }

    /// Compute `a^2` (mod l)
    #[inline(never)]
    #[allow(dead_code)] // XXX we don't expose square() via the Scalar API
    pub(crate) fn square(&self) -> Scalar52 {
        let aa = Scalar52::montgomery_reduce(&Scalar52::square_internal(self));
        Scalar52::montgomery_reduce(&Scalar52::mul_internal(&aa, &constants::RR))
    }
}
