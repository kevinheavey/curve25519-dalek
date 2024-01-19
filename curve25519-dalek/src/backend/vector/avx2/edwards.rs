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

//! Parallel Edwards Arithmetic for Curve25519.
//!
//! This module currently has two point types:
//!
//! * `ExtendedPoint`: a point stored in vector-friendly format, with
//! vectorized doubling and addition;
//!
//! * `CachedPoint`: used for readdition.
//!
//! Details on the formulas can be found in the documentation for the
//! parent `avx2` module.
//!
//! This API is designed to be safe: vectorized points can only be
//! created from serial points (which do validation on decompression),
//! and operations on valid points return valid points, so invalid
//! point states should be unrepresentable.
//!
//! This design goal is met, with one exception: the `Neg`
//! implementation for the `CachedPoint` performs a lazy negation, so
//! that subtraction can be efficiently implemented as a negation and
//! an addition.  Repeatedly negating a `CachedPoint` will cause its
//! coefficients to grow and eventually overflow.  Repeatedly negating
//! a point should not be necessary anyways.

#![allow(non_snake_case)]


use super::field::{FieldElement2625x4};

/// A point on Curve25519, using parallel Edwards formulas for curve
/// operations.
///
/// # Invariant
///
/// The coefficients of an `ExtendedPoint` are bounded with
/// \\( b < 0.007 \\).
#[derive(Copy, Clone, Debug)]
pub(crate) struct ExtendedPoint(pub(super) FieldElement2625x4);


/// A cached point with some precomputed variables used for readdition.
///
/// # Warning
///
/// It is not safe to negate this point more than once.
///
/// # Invariant
///
/// As long as the `CachedPoint` is not repeatedly negated, its
/// coefficients will be bounded with \\( b < 1.0 \\).
#[derive(Copy, Clone, Debug)]
pub(crate) struct CachedPoint(pub(super) FieldElement2625x4);
