// -*- mode: rust; -*-
//
// This file is part of curve25519-dalek.
// Copyright (c) 2016-2021 isis agora lovecruft
// Copyright (c) 2016-2019 Henry de Valence
// See LICENSE for licensing information.
//
// Authors:
// - Isis Agora Lovecruft <isis@patternsinthevoid.net>
// - Henry de Valence <hdevalence@hdevalence.ca>

//! Field arithmetic modulo \\(p = 2\^{255} - 19\\).
//!
//! The `curve25519_dalek::field` module provides a type alias
//! `curve25519_dalek::field::FieldElement` to a field element type
//! defined in the `backend` module; either `FieldElement51` or
//! `FieldElement2625`.
//!
//! Field operations defined in terms of machine
//! operations, such as field multiplication or squaring, are defined in
//! the backend implementation.
//!
//! Field operations defined in terms of other field operations, such as
//! field inversion or square roots, are defined here.

#![allow(unused_qualifications)]

use core::cmp::{Eq, PartialEq};

use cfg_if::cfg_if;

use subtle::Choice;
use subtle::ConstantTimeEq;

use crate::backend;

cfg_if! {
    if #[cfg(curve25519_dalek_backend = "fiat")] {
        /// A `FieldElement` represents an element of the field
        /// \\( \mathbb Z / (2\^{255} - 19)\\).
        ///
        /// The `FieldElement` type is an alias for one of the platform-specific
        /// implementations.
        ///
        /// Using formally-verified field arithmetic from fiat-crypto.
        #[cfg(curve25519_dalek_bits = "32")]
        pub(crate) type FieldElement = backend::serial::fiat_u32::field::FieldElement2625;

        /// A `FieldElement` represents an element of the field
        /// \\( \mathbb Z / (2\^{255} - 19)\\).
        ///
        /// The `FieldElement` type is an alias for one of the platform-specific
        /// implementations.
        ///
        /// Using formally-verified field arithmetic from fiat-crypto.
        #[cfg(curve25519_dalek_bits = "64")]
        pub(crate) type FieldElement = backend::serial::fiat_u64::field::FieldElement51;
    } else if #[cfg(curve25519_dalek_bits = "64")] {
        /// A `FieldElement` represents an element of the field
        /// \\( \mathbb Z / (2\^{255} - 19)\\).
        ///
        /// The `FieldElement` type is an alias for one of the platform-specific
        /// implementations.
        pub(crate) type FieldElement = backend::serial::u64::field::FieldElement51;
    } else {
        /// A `FieldElement` represents an element of the field
        /// \\( \mathbb Z / (2\^{255} - 19)\\).
        ///
        /// The `FieldElement` type is an alias for one of the platform-specific
        /// implementations.
        pub(crate) type FieldElement = backend::serial::u32::field::FieldElement2625;
    }
}

impl Eq for FieldElement {}

impl PartialEq for FieldElement {
    fn eq(&self, other: &FieldElement) -> bool {
        self.ct_eq(other).into()
    }
}

impl ConstantTimeEq for FieldElement {
    /// Test equality between two `FieldElement`s.  Since the
    /// internal representation is not canonical, the field elements
    /// are normalized to wire format before comparison.
    fn ct_eq(&self, other: &FieldElement) -> Choice {
        self.as_bytes().ct_eq(&other.as_bytes())
    }
}

impl FieldElement {
    /// Compute (self^(2^250-1), self^11), used as a helper function
    /// within invert() and pow22523().
    #[rustfmt::skip] // keep alignment of explanatory comments
    fn pow22501(&self) -> (FieldElement, FieldElement) {
        // Instead of managing which temporary variables are used
        // for what, we define as many as we need and leave stack
        // allocation to the compiler
        //
        // Each temporary variable t_i is of the form (self)^e_i.
        // Squaring t_i corresponds to multiplying e_i by 2,
        // so the pow2k function shifts e_i left by k places.
        // Multiplying t_i and t_j corresponds to adding e_i + e_j.
        //
        // Temporary t_i                      Nonzero bits of e_i
        //
        let t0  = self.square();           // 1         e_0 = 2^1
        let t1  = t0.square().square();    // 3         e_1 = 2^3
        let t2  = self * &t1;              // 3,0       e_2 = 2^3 + 2^0
        let t3  = &t0 * &t2;               // 3,1,0
        let t4  = t3.square();             // 4,2,1
        let t5  = &t2 * &t4;               // 4,3,2,1,0
        let t6  = t5.pow2k(5);             // 9,8,7,6,5
        let t7  = &t6 * &t5;               // 9,8,7,6,5,4,3,2,1,0
        let t8  = t7.pow2k(10);            // 19..10
        let t9  = &t8 * &t7;               // 19..0
        let t10 = t9.pow2k(20);            // 39..20
        let t11 = &t10 * &t9;              // 39..0
        let t12 = t11.pow2k(10);           // 49..10
        let t13 = &t12 * &t7;              // 49..0
        let t14 = t13.pow2k(50);           // 99..50
        let t15 = &t14 * &t13;             // 99..0
        let t16 = t15.pow2k(100);          // 199..100
        let t17 = &t16 * &t15;             // 199..0
        let t18 = t17.pow2k(50);           // 249..50
        let t19 = &t18 * &t13;             // 249..0

        (t19, t3)
    }



    /// Raise this field element to the power (p-5)/8 = 2^252 -3.
    #[rustfmt::skip] // keep alignment of explanatory comments
    #[allow(clippy::let_and_return)]
    fn pow_p58(&self) -> FieldElement {
        // The bits of (p-5)/8 are 101111.....11.
        //
        //                                 nonzero bits of exponent
        let (t19, _) = self.pow22501();    // 249..0
        let t20 = t19.pow2k(2);            // 251..2
        let t21 = self * &t20;             // 251..2,0

        t21
    }

    pub(crate) fn was_nonzero_square(u: &FieldElement, v: &FieldElement) -> Choice {
        // Using the same trick as in ed25519 decoding, we merge the
        // inversion, the square root, and the square test as follows.
        //
        // To compute sqrt(α), we can compute β = α^((p+3)/8).
        // Then β^2 = ±α, so multiplying β by sqrt(-1) if necessary
        // gives sqrt(α).
        //
        // To compute 1/sqrt(α), we observe that
        //    1/β = α^(p-1 - (p+3)/8) = α^((7p-11)/8)
        //                            = α^3 * (α^7)^((p-5)/8).
        //
        // We can therefore compute sqrt(u/v) = sqrt(u)/sqrt(v)
        // by first computing
        //    r = u^((p+3)/8) v^(p-1-(p+3)/8)
        //      = u u^((p-5)/8) v^3 (v^7)^((p-5)/8)
        //      = (uv^3) (uv^7)^((p-5)/8).
        //
        // If v is nonzero and u/v is square, then r^2 = ±u/v,
        //                                     so vr^2 = ±u.
        // If vr^2 =  u, then sqrt(u/v) = r.
        // If vr^2 = -u, then sqrt(u/v) = r*sqrt(-1).
        //
        // If v is zero, r is also zero.

        let v3 = &v.square() * v;
        let v7 = &v3.square() * v;
        let r = &(u * &v3) * &(u * &v7).pow_p58();
        let check = v * &r.square();
        let correct_sign_sqrt = check.ct_eq(u);
        let flipped_sign_sqrt = check.ct_eq(&(-u));
        correct_sign_sqrt | flipped_sign_sqrt
    }
}
