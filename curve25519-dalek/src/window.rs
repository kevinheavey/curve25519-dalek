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

//! Code for fixed- and sliding-window functionality

#![allow(non_snake_case)]

use core::fmt::Debug;

use cfg_if::cfg_if;

use crate::backend::serial::curve_models::AffineNielsPoint;
use crate::backend::serial::curve_models::ProjectiveNielsPoint;
use crate::edwards::EdwardsPoint;

macro_rules! impl_lookup_table {
    (Name = $name:ident, Size = $size:expr, SizeNeg = $neg:expr, SizeRange = $range:expr, ConversionRange = $conv_range:expr) => {
        /// A lookup table of precomputed multiples of a point \\(P\\), used to
        /// compute \\( xP \\) for \\( -8 \leq x \leq 8 \\).
        ///
        /// The computation of \\( xP \\) is done in constant time by the `select` function.
        ///
        /// Since `LookupTable` does not implement `Index`, it's more difficult
        /// to accidentally use the table directly.  Unfortunately the table is
        /// only `pub(crate)` so that we can write hardcoded constants, so it's
        /// still technically possible.  It would be nice to prevent direct
        /// access to the table.
        #[derive(Copy, Clone)]
        pub(crate) struct $name<T>(pub(crate) [T; $size]);
    };
} // End macro_rules! impl_lookup_table

// The first one has to be named "LookupTable" because it's used as a constructor for consts.
// This is radix-16
impl_lookup_table! {
    Name = LookupTable,
    Size = 8,
    SizeNeg = -8,
    SizeRange = 1..9,
    ConversionRange = 0..7
}

// The rest only get used to make basepoint tables
cfg_if! {
    if #[cfg(feature = "precomputed-tables")] {
        // radix-32
        impl_lookup_table! {
            Name = LookupTableRadix32,
            Size = 16,
            SizeNeg = -16,
            SizeRange = 1..17,
            ConversionRange = 0..15
        }
        // radix-64
        impl_lookup_table! {
            Name = LookupTableRadix64,
            Size = 32,
            SizeNeg = -32,
            SizeRange = 1..33,
            ConversionRange = 0..31
        }
        // radix-128
        impl_lookup_table! {
            Name = LookupTableRadix128,
            Size = 64,
            SizeNeg = -64,
            SizeRange = 1..65,
            ConversionRange = 0..63
        }
        // radix-256
        impl_lookup_table! {
            Name = LookupTableRadix256,
            Size = 128,
            SizeNeg = -128,
            SizeRange = 1..129,
            ConversionRange = 0..127
        }

        // For homogeneity we then alias it to "LookupTableRadix16".
        pub(crate) type LookupTableRadix16<T> = LookupTable<T>;
    }
}

/// Holds odd multiples 1A, 3A, ..., 15A of a point A.
#[derive(Copy, Clone)]
pub(crate) struct NafLookupTable5<T>(pub(crate) [T; 8]);

impl<T: Debug> Debug for NafLookupTable5<T> {
    fn fmt(&self, f: &mut ::core::fmt::Formatter<'_>) -> ::core::fmt::Result {
        write!(f, "NafLookupTable5({:?})", self.0)
    }
}

impl<'a> From<&'a EdwardsPoint> for NafLookupTable5<ProjectiveNielsPoint> {
    fn from(A: &'a EdwardsPoint) -> Self {
        let mut Ai = [A.as_projective_niels(); 8];
        let A2 = A.double();
        for i in 0..7 {
            Ai[i + 1] = (&A2 + &Ai[i]).as_extended().as_projective_niels();
        }
        // Now Ai = [A, 3A, 5A, 7A, 9A, 11A, 13A, 15A]
        NafLookupTable5(Ai)
    }
}

impl<'a> From<&'a EdwardsPoint> for NafLookupTable5<AffineNielsPoint> {
    fn from(A: &'a EdwardsPoint) -> Self {
        let mut Ai = [A.as_affine_niels(); 8];
        let A2 = A.double();
        for i in 0..7 {
            Ai[i + 1] = (&A2 + &Ai[i]).as_extended().as_affine_niels();
        }
        // Now Ai = [A, 3A, 5A, 7A, 9A, 11A, 13A, 15A]
        NafLookupTable5(Ai)
    }
}

/// Holds stuff up to 8. The only time we use tables this big is for precomputed basepoint tables
/// and multiscalar multiplication (which requires alloc).
#[cfg(any(feature = "precomputed-tables", feature = "alloc"))]
#[derive(Copy, Clone)]
pub(crate) struct NafLookupTable8<T>(pub(crate) [T; 64]);


#[cfg(any(feature = "precomputed-tables", feature = "alloc"))]
impl<T: Debug> Debug for NafLookupTable8<T> {
    fn fmt(&self, f: &mut ::core::fmt::Formatter<'_>) -> ::core::fmt::Result {
        writeln!(f, "NafLookupTable8([")?;
        for i in 0..64 {
            writeln!(f, "\t{:?},", &self.0[i])?;
        }
        write!(f, "])")
    }
}

#[cfg(any(feature = "precomputed-tables", feature = "alloc"))]
impl<'a> From<&'a EdwardsPoint> for NafLookupTable8<ProjectiveNielsPoint> {
    fn from(A: &'a EdwardsPoint) -> Self {
        let mut Ai = [A.as_projective_niels(); 64];
        let A2 = A.double();
        for i in 0..63 {
            Ai[i + 1] = (&A2 + &Ai[i]).as_extended().as_projective_niels();
        }
        // Now Ai = [A, 3A, 5A, 7A, 9A, 11A, 13A, 15A, ..., 127A]
        NafLookupTable8(Ai)
    }
}

#[cfg(any(feature = "precomputed-tables", feature = "alloc"))]
impl<'a> From<&'a EdwardsPoint> for NafLookupTable8<AffineNielsPoint> {
    fn from(A: &'a EdwardsPoint) -> Self {
        let mut Ai = [A.as_affine_niels(); 64];
        let A2 = A.double();
        for i in 0..63 {
            Ai[i + 1] = (&A2 + &Ai[i]).as_extended().as_affine_niels();
        }
        // Now Ai = [A, 3A, 5A, 7A, 9A, 11A, 13A, 15A, ..., 127A]
        NafLookupTable8(Ai)
    }
}
