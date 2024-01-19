// -*- mode: rust; -*-
//
// This file is part of curve25519-dalek.
// Copyright (c) 2016-2021 isis lovecruft
// Copyright (c) 2016-2019 Henry de Valence
// See LICENSE for licensing information.
//
// Authors:
// - isis agora lovecruft <isis@patternsinthevoid.net>
// - Henry de Valence <hdevalence@hdevalence.ca>

//! Scalar multiplication on the Montgomery form of Curve25519.
//!
//! To avoid notational confusion with the Edwards code, we use
//! variables \\( u, v \\) for the Montgomery curve, so that “Montgomery
//! \\(u\\)” here corresponds to “Montgomery \\(x\\)” elsewhere.
//!
//! Montgomery arithmetic works not on the curve itself, but on the
//! \\(u\\)-line, which discards sign information and unifies the curve
//! and its quadratic twist.  See [_Montgomery curves and their
//! arithmetic_][costello-smith] by Costello and Smith for more details.
//!
//! The `MontgomeryPoint` struct contains the affine \\(u\\)-coordinate
//! \\(u\_0(P)\\) of a point \\(P\\) on either the curve or the twist.
//! Here the map \\(u\_0 : \mathcal M \rightarrow \mathbb F\_p \\) is
//! defined by \\(u\_0((u,v)) = u\\); \\(u\_0(\mathcal O) = 0\\).  See
//! section 5.4 of Costello-Smith for more details.
//!
//! # Scalar Multiplication
//!
//! Scalar multiplication on `MontgomeryPoint`s is provided by the `*`
//! operator, which implements the Montgomery ladder.
//!
//! # Edwards Conversion
//!
//! The \\(2\\)-to-\\(1\\) map from the Edwards model to the Montgomery
//! \\(u\\)-line is provided by `EdwardsPoint::to_montgomery()`.
//!
//! To lift a `MontgomeryPoint` to an `EdwardsPoint`, use
//! `MontgomeryPoint::to_edwards()`, which takes a sign parameter.
//! This function rejects `MontgomeryPoints` which correspond to points
//! on the twist.
//!
//! [costello-smith]: https://eprint.iacr.org/2017/212.pdf

// We allow non snake_case names because coordinates in projective space are
// traditionally denoted by the capitalisation of their respective
// counterparts in affine space.  Yeah, you heard me, rustc, I'm gonna have my
// affine and projective cakes and eat both of them too.
#![allow(non_snake_case)]

use core::{
    hash::{Hash, Hasher},
    ops::{Mul, MulAssign},
};

use crate::constants::{APLUS2_OVER_FOUR, MONTGOMERY_A, MONTGOMERY_A_NEG};
use crate::field::FieldElement;
use crate::scalar::{Scalar};

use crate::traits::Identity;

use subtle::Choice;
use subtle::ConstantTimeEq;
use subtle::{ConditionallyNegatable, ConditionallySelectable};

#[cfg(feature = "zeroize")]
use zeroize::Zeroize;

/// Holds the \\(u\\)-coordinate of a point on the Montgomery form of
/// Curve25519 or its twist.
#[derive(Copy, Clone, Debug, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub(crate) struct MontgomeryPoint(pub [u8; 32]);

/// Equality of `MontgomeryPoint`s is defined mod p.
impl ConstantTimeEq for MontgomeryPoint {
    fn ct_eq(&self, other: &MontgomeryPoint) -> Choice {
        let self_fe = FieldElement::from_bytes(&self.0);
        let other_fe = FieldElement::from_bytes(&other.0);

        self_fe.ct_eq(&other_fe)
    }
}

impl PartialEq for MontgomeryPoint {
    fn eq(&self, other: &MontgomeryPoint) -> bool {
        self.ct_eq(other).into()
    }
}

impl Eq for MontgomeryPoint {}

// Equal MontgomeryPoints must hash to the same value. So we have to get them into a canonical
// encoding first
impl Hash for MontgomeryPoint {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // Do a round trip through a `FieldElement`. `as_bytes` is guaranteed to give a canonical
        // 32-byte encoding
        let canonical_bytes = FieldElement::from_bytes(&self.0).as_bytes();
        canonical_bytes.hash(state);
    }
}

impl Identity for MontgomeryPoint {
    /// Return the group identity element, which has order 4.
    fn identity() -> MontgomeryPoint {
        MontgomeryPoint([0u8; 32])
    }
}

#[cfg(feature = "zeroize")]
impl Zeroize for MontgomeryPoint {
    fn zeroize(&mut self) {
        self.0.zeroize();
    }
}

impl MontgomeryPoint {
    /// Given `self` \\( = u\_0(P) \\), and a big-endian bit representation of an integer
    /// \\(n\\), return \\( u\_0(\[n\]P) \\). This is constant time in the length of `bits`.
    ///
    /// **NOTE:** You probably do not want to use this function. Almost every protocol built on
    /// Curve25519 uses _clamped multiplication_, explained
    /// [here](https://neilmadden.blog/2020/05/28/whats-the-curve25519-clamping-all-about/).
    /// When in doubt, use [`Self::mul_clamped`].
    pub(crate) fn mul_bits_be(&self, bits: impl Iterator<Item = bool>) -> MontgomeryPoint {
        // Algorithm 8 of Costello-Smith 2017
        let affine_u = FieldElement::from_bytes(&self.0);
        let mut x0 = ProjectivePoint::identity();
        let mut x1 = ProjectivePoint {
            U: affine_u,
            W: FieldElement::ONE,
        };

        // Go through the bits from most to least significant, using a sliding window of 2
        let mut prev_bit = false;
        for cur_bit in bits {
            let choice: u8 = (prev_bit ^ cur_bit) as u8;

            debug_assert!(choice == 0 || choice == 1);

            ProjectivePoint::conditional_swap(&mut x0, &mut x1, choice.into());
            differential_add_and_double(&mut x0, &mut x1, &affine_u);

            prev_bit = cur_bit;
        }
        // The final value of prev_bit above is scalar.bits()[0], i.e., the LSB of scalar
        ProjectivePoint::conditional_swap(&mut x0, &mut x1, Choice::from(prev_bit as u8));
        // Don't leave the bit in the stack
        #[cfg(feature = "zeroize")]
        prev_bit.zeroize();

        x0.as_affine()
    }
}

/// Perform the Elligator2 mapping to a Montgomery point.
///
/// See <https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-10#section-6.7.1>
//
// TODO Determine how much of the hash-to-group API should be exposed after the CFRG
//      draft gets into a more polished/accepted state.
#[allow(unused)]
pub(crate) fn elligator_encode(r_0: &FieldElement) -> MontgomeryPoint {
    let one = FieldElement::ONE;
    let d_1 = &one + &r_0.square2(); /* 2r^2 */

    let d = &MONTGOMERY_A_NEG * &(d_1.invert()); /* A/(1+2r^2) */

    let d_sq = &d.square();
    let au = &MONTGOMERY_A * &d;

    let inner = &(d_sq + &au) + &one;
    let eps = &d * &inner; /* eps = d^3 + Ad^2 + d */

    let (eps_is_sq, _eps) = FieldElement::sqrt_ratio_i(&eps, &one);

    let zero = FieldElement::ZERO;
    let Atemp = FieldElement::conditional_select(&MONTGOMERY_A, &zero, eps_is_sq); /* 0, or A if nonsquare*/
    let mut u = &d + &Atemp; /* d, or d+A if nonsquare */
    u.conditional_negate(!eps_is_sq); /* d, or -d-A if nonsquare */

    MontgomeryPoint(u.as_bytes())
}

/// A `ProjectivePoint` holds a point on the projective line
/// \\( \mathbb P(\mathbb F\_p) \\), which we identify with the Kummer
/// line of the Montgomery curve.
#[derive(Copy, Clone, Debug)]
struct ProjectivePoint {
    pub U: FieldElement,
    pub W: FieldElement,
}

impl Identity for ProjectivePoint {
    fn identity() -> ProjectivePoint {
        ProjectivePoint {
            U: FieldElement::ONE,
            W: FieldElement::ZERO,
        }
    }
}

impl Default for ProjectivePoint {
    fn default() -> ProjectivePoint {
        ProjectivePoint::identity()
    }
}

impl ConditionallySelectable for ProjectivePoint {
    fn conditional_select(
        a: &ProjectivePoint,
        b: &ProjectivePoint,
        choice: Choice,
    ) -> ProjectivePoint {
        ProjectivePoint {
            U: FieldElement::conditional_select(&a.U, &b.U, choice),
            W: FieldElement::conditional_select(&a.W, &b.W, choice),
        }
    }
}

impl ProjectivePoint {
    /// Dehomogenize this point to affine coordinates.
    ///
    /// # Return
    ///
    /// * \\( u = U / W \\) if \\( W \neq 0 \\);
    /// * \\( 0 \\) if \\( W \eq 0 \\);
    pub(crate) fn as_affine(&self) -> MontgomeryPoint {
        let u = &self.U * &self.W.invert();
        MontgomeryPoint(u.as_bytes())
    }
}

/// Perform the double-and-add step of the Montgomery ladder.
///
/// Given projective points
/// \\( (U\_P : W\_P) = u(P) \\),
/// \\( (U\_Q : W\_Q) = u(Q) \\),
/// and the affine difference
/// \\(      u\_{P-Q} = u(P-Q) \\), set
/// $$
///     (U\_P : W\_P) \gets u(\[2\]P)
/// $$
/// and
/// $$
///     (U\_Q : W\_Q) \gets u(P + Q).
/// $$
#[rustfmt::skip] // keep alignment of explanatory comments
fn differential_add_and_double(
    P: &mut ProjectivePoint,
    Q: &mut ProjectivePoint,
    affine_PmQ: &FieldElement,
) {
    let t0 = &P.U + &P.W;
    let t1 = &P.U - &P.W;
    let t2 = &Q.U + &Q.W;
    let t3 = &Q.U - &Q.W;

    let t4 = t0.square();   // (U_P + W_P)^2 = U_P^2 + 2 U_P W_P + W_P^2
    let t5 = t1.square();   // (U_P - W_P)^2 = U_P^2 - 2 U_P W_P + W_P^2

    let t6 = &t4 - &t5;     // 4 U_P W_P

    let t7 = &t0 * &t3;     // (U_P + W_P) (U_Q - W_Q) = U_P U_Q + W_P U_Q - U_P W_Q - W_P W_Q
    let t8 = &t1 * &t2;     // (U_P - W_P) (U_Q + W_Q) = U_P U_Q - W_P U_Q + U_P W_Q - W_P W_Q

    let t9  = &t7 + &t8;    // 2 (U_P U_Q - W_P W_Q)
    let t10 = &t7 - &t8;    // 2 (W_P U_Q - U_P W_Q)

    let t11 =  t9.square(); // 4 (U_P U_Q - W_P W_Q)^2
    let t12 = t10.square(); // 4 (W_P U_Q - U_P W_Q)^2

    let t13 = &APLUS2_OVER_FOUR * &t6; // (A + 2) U_P U_Q

    let t14 = &t4 * &t5;    // ((U_P + W_P)(U_P - W_P))^2 = (U_P^2 - W_P^2)^2
    let t15 = &t13 + &t5;   // (U_P - W_P)^2 + (A + 2) U_P W_P

    let t16 = &t6 * &t15;   // 4 (U_P W_P) ((U_P - W_P)^2 + (A + 2) U_P W_P)

    let t17 = affine_PmQ * &t12; // U_D * 4 (W_P U_Q - U_P W_Q)^2
    let t18 = t11;               // W_D * 4 (U_P U_Q - W_P W_Q)^2

    P.U = t14;  // U_{P'} = (U_P + W_P)^2 (U_P - W_P)^2
    P.W = t16;  // W_{P'} = (4 U_P W_P) ((U_P - W_P)^2 + ((A + 2)/4) 4 U_P W_P)
    Q.U = t18;  // U_{Q'} = W_D * 4 (U_P U_Q - W_P W_Q)^2
    Q.W = t17;  // W_{Q'} = U_D * 4 (W_P U_Q - U_P W_Q)^2
}

define_mul_assign_variants!(LHS = MontgomeryPoint, RHS = Scalar);

define_mul_variants!(
    LHS = MontgomeryPoint,
    RHS = Scalar,
    Output = MontgomeryPoint
);
define_mul_variants!(
    LHS = Scalar,
    RHS = MontgomeryPoint,
    Output = MontgomeryPoint
);

/// Multiply this `MontgomeryPoint` by a `Scalar`.
impl Mul<&Scalar> for &MontgomeryPoint {
    type Output = MontgomeryPoint;

    /// Given `self` \\( = u\_0(P) \\), and a `Scalar` \\(n\\), return \\( u\_0(\[n\]P) \\)
    fn mul(self, scalar: &Scalar) -> MontgomeryPoint {
        // We multiply by the integer representation of the given Scalar. By scalar invariant #1,
        // the MSB is 0, so we can skip it.
        self.mul_bits_be(scalar.bits_le().rev().skip(1))
    }
}

impl MulAssign<&Scalar> for MontgomeryPoint {
    fn mul_assign(&mut self, scalar: &Scalar) {
        *self = (self as &MontgomeryPoint) * scalar;
    }
}

impl Mul<&MontgomeryPoint> for &Scalar {
    type Output = MontgomeryPoint;

    fn mul(self, point: &MontgomeryPoint) -> MontgomeryPoint {
        point * self
    }
}
