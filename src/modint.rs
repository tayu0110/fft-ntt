use num::{
    traits::{
        ops::overflowing::{OverflowingAdd, OverflowingSub},
        WrappingAdd,
    },
    Integer, One, PrimInt, Zero,
};
use std::marker;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub, SubAssign};

pub trait Modulo<T = i64>: Clone + marker::Copy + PartialEq + Eq {
    fn modulo() -> T;
    fn primitive_root() -> T {
        unimplemented!()
    }
    fn montgomery_constant_mask() -> T {
        unimplemented!()
    }
    fn montgomery_constant_r() -> T {
        unimplemented!()
    }
    fn montgomery_constant_r_inv() -> T {
        unimplemented!()
    }
    fn montgomery_constant_r_pow2() -> T {
        unimplemented!()
    }
    fn montgomery_constant_modulo_inv() -> T {
        unimplemented!()
    }
}

#[derive(Clone, marker::Copy, PartialEq, Eq)]
pub enum Mod998244353<T = i64> {
    PhantomData(std::marker::PhantomData<T>),
}
impl Modulo for Mod998244353 {
    #[inline]
    fn modulo() -> i64 {
        998_244_353i64
    }
    #[inline]
    // R - 1 = 2^63 - 1
    fn montgomery_constant_mask() -> i64 {
        0x7FFFFFFFFFFFFFFF
    }
    #[inline]
    fn primitive_root() -> i64 {
        3i64
    }
    #[inline]
    // R = 2^63 mod 998244353
    fn montgomery_constant_r() -> i64 {
        466025955
    }
    #[inline]
    // R2 = 2^126 mod 998244353
    fn montgomery_constant_r_pow2() -> i64 {
        74890016
    }
    #[inline]
    // R^{-1} = (2^63 mod 998244353)^{-1} mod 998244353
    fn montgomery_constant_r_inv() -> i64 {
        890394177
    }
    #[inline]
    // modulo * modulo_inv = -1 mod R
    fn montgomery_constant_modulo_inv() -> i64 {
        8226880251553120255
    }
}
impl Modulo<u32> for Mod998244353<u32> {
    #[inline]
    fn modulo() -> u32 {
        998_244_353u32
    }
    #[inline]
    // R - 1 = 2^32 - 1
    fn montgomery_constant_mask() -> u32 {
        !0
    }
    #[inline]
    // modulo * modulo_inv = -1 mod R
    fn montgomery_constant_modulo_inv() -> u32 {
        998244351
    }
    #[inline]
    // R = 2^32 mod 998244353
    fn montgomery_constant_r() -> u32 {
        301989884
    }
    #[inline]
    // R^{-1} = (2^32 mod 998244353)^{-1} mod 998244353
    fn montgomery_constant_r_inv() -> u32 {
        232013824
    }
    #[inline]
    // R2 = 2^64 mod 998244353
    fn montgomery_constant_r_pow2() -> u32 {
        932051910
    }
    #[inline]
    fn primitive_root() -> u32 {
        3
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// MontgomeryModint
////////////////////////////////////////////////////////////////////////////////////////////////////////////
pub trait MontgomeryMultiplication<M: Modulo<T>, T> {
    // t <- MR(T) = (T + (TN' mod R) * N) / R
    //  if t >= N then return t - N else return t
    //      T := a (0 <= T < NR)
    //      N := modulo()
    //      N':= montgomery_constant_modulo_inv()
    //      R := montgomery_constant_r()
    fn montgomery_reduction(self) -> Self;
    fn multiplication(self, rhs: Self) -> Self;
}

impl<M: Modulo<u32>> MontgomeryMultiplication<M, u32> for u32 {
    #[inline]
    fn montgomery_reduction(self) -> u32 {
        let t = ((self as u64).wrapping_add(
            (self.wrapping_mul(M::montgomery_constant_modulo_inv()) as u64)
                .wrapping_mul(M::modulo() as u64),
        ) >> 32) as u32;
        if t >= M::modulo() {
            t - M::modulo()
        } else {
            t
        }
    }
    #[inline]
    fn multiplication(self, rhs: Self) -> Self {
        let a = self as u64 * rhs as u64;
        let t = (a.wrapping_add(
            ((a as u32).wrapping_mul(M::montgomery_constant_modulo_inv()) as u64)
                .wrapping_mul(M::modulo() as u64),
        ) >> 32) as u32;
        if t >= M::modulo() {
            t - M::modulo()
        } else {
            t
        }
    }
}

impl<M: Modulo> MontgomeryMultiplication<M, i64> for i64 {
    #[inline]
    fn montgomery_reduction(self) -> i64 {
        let t = ((self as i128).wrapping_add(
            ((self as i128).wrapping_mul(M::montgomery_constant_modulo_inv() as i128)
                & M::montgomery_constant_mask() as i128)
                .wrapping_mul(M::modulo() as i128),
        ) >> 63) as i64;
        if t >= M::modulo() {
            t - M::modulo()
        } else {
            t
        }
    }
    #[inline]
    fn multiplication(self, rhs: Self) -> Self {
        let a = self as i128 * rhs as i128;
        let t = (a.wrapping_add(
            (a.wrapping_mul(M::montgomery_constant_modulo_inv() as i128)
                & M::montgomery_constant_mask() as i128)
                .wrapping_mul(M::modulo() as i128),
        ) >> 63) as i64;
        if t >= M::modulo() {
            t - M::modulo()
        } else {
            t
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub struct MontgomeryModint<M: Modulo<T>, T = i64>
where
    T: Integer + PrimInt + MontgomeryMultiplication<M, T>,
{
    val: T,
    _phantom: std::marker::PhantomData<M>,
}

impl<
        M: Modulo<T>,
        T: Integer
            + PrimInt
            + MontgomeryMultiplication<M, T>
            + OverflowingAdd
            + OverflowingSub
            + WrappingAdd,
    > MontgomeryModint<M, T>
{
    #[inline]
    pub fn new(val: T) -> Self {
        let val = val.multiplication(M::montgomery_constant_r_pow2());
        Self {
            val,
            _phantom: std::marker::PhantomData,
        }
    }

    #[inline]
    pub fn val(&self) -> T {
        self.val.montgomery_reduction()
    }

    #[inline]
    pub fn val_montgomery_expression(&self) -> T {
        self.val
    }

    #[inline]
    pub fn one() -> Self {
        Self {
            val: M::montgomery_constant_r(),
            _phantom: std::marker::PhantomData,
        }
    }

    #[inline]
    pub fn zero() -> Self {
        Self {
            val: T::zero(),
            _phantom: std::marker::PhantomData,
        }
    }

    pub fn pow(&self, mut n: T) -> Self {
        let mut val = self.val;
        let mut res = if (n & T::one()) != T::zero() {
            val
        } else {
            M::montgomery_constant_r()
        };
        n = n.signed_shr(1);
        while n != T::zero() {
            val = val.multiplication(val);
            if n & T::one() != T::zero() {
                res = res.multiplication(val);
            }
            n = n.signed_shr(1);
        }
        Self {
            val: res,
            _phantom: std::marker::PhantomData,
        }
    }

    #[inline]
    pub fn nth_root(n: T) -> Self {
        debug_assert!(n == T::one().signed_shl(n.trailing_zeros()));
        debug_assert!(M::modulo() - T::one() + (M::modulo() - T::one()) / n >= T::zero());

        MontgomeryModint::<M, T>::new(M::primitive_root())
            .pow(M::modulo() - T::one() + (M::modulo() - T::one()) / n)
    }

    #[inline]
    pub fn inv(&self) -> Self {
        self.pow(M::modulo() - T::one() - T::one())
    }
}

impl<
        M: Modulo<T>,
        T: Integer
            + PrimInt
            + MontgomeryMultiplication<M, T>
            + OverflowingAdd
            + OverflowingSub
            + WrappingAdd,
    > One for MontgomeryModint<M, T>
{
    #[inline]
    fn one() -> Self {
        Self::one()
    }
}

impl<
        M: Modulo<T>,
        T: Integer
            + PrimInt
            + MontgomeryMultiplication<M, T>
            + OverflowingAdd
            + OverflowingSub
            + WrappingAdd,
    > Zero for MontgomeryModint<M, T>
{
    #[inline]
    fn zero() -> Self {
        Self::zero()
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.val == Self::zero().val
    }
}

impl<
        M: Modulo<T>,
        T: Integer
            + PrimInt
            + MontgomeryMultiplication<M, T>
            + OverflowingAdd
            + OverflowingSub
            + WrappingAdd,
    > Add for MontgomeryModint<M, T>
{
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        let (t, fa) = self.val.overflowing_add(&rhs.val);
        let (u, fs) = t.overflowing_sub(&M::modulo());
        Self {
            val: if fa || !fs { u } else { t },
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<
        M: Modulo<T>,
        T: Integer
            + PrimInt
            + MontgomeryMultiplication<M, T>
            + OverflowingAdd
            + OverflowingSub
            + WrappingAdd,
    > Sub for MontgomeryModint<M, T>
{
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        let (val, f) = self.val.overflowing_sub(&rhs.val);
        Self {
            val: if f || val < T::zero() {
                val.wrapping_add(&M::modulo())
            } else {
                val
            },
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<
        M: Modulo<T>,
        T: Integer
            + PrimInt
            + MontgomeryMultiplication<M, T>
            + OverflowingAdd
            + OverflowingSub
            + WrappingAdd,
    > Mul for MontgomeryModint<M, T>
{
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        Self {
            val: self.val.multiplication(rhs.val),
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<
        M: Modulo<T>,
        T: Integer
            + PrimInt
            + MontgomeryMultiplication<M, T>
            + OverflowingAdd
            + OverflowingSub
            + WrappingAdd,
    > Div for MontgomeryModint<M, T>
{
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        self * rhs.inv()
    }
}

impl<
        M: Modulo<T>,
        T: Integer
            + PrimInt
            + MontgomeryMultiplication<M, T>
            + OverflowingAdd
            + OverflowingSub
            + WrappingAdd,
    > AddAssign for MontgomeryModint<M, T>
{
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl<
        M: Modulo<T>,
        T: Integer
            + PrimInt
            + MontgomeryMultiplication<M, T>
            + OverflowingAdd
            + OverflowingSub
            + WrappingAdd,
    > SubAssign for MontgomeryModint<M, T>
{
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl<
        M: Modulo<T>,
        T: Integer
            + PrimInt
            + MontgomeryMultiplication<M, T>
            + OverflowingAdd
            + OverflowingSub
            + WrappingAdd,
    > MulAssign for MontgomeryModint<M, T>
{
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl<
        M: Modulo<T>,
        T: Integer
            + PrimInt
            + MontgomeryMultiplication<M, T>
            + OverflowingAdd
            + OverflowingSub
            + WrappingAdd,
    > DivAssign for MontgomeryModint<M, T>
{
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs;
    }
}

impl<
        M: Modulo<T>,
        T: Integer
            + PrimInt
            + MontgomeryMultiplication<M, T>
            + OverflowingAdd
            + OverflowingSub
            + WrappingAdd
            + std::fmt::Display,
    > std::fmt::Debug for MontgomeryModint<M, T>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.val())
    }
}

impl<
        M: Modulo<T>,
        T: Integer
            + PrimInt
            + MontgomeryMultiplication<M, T>
            + OverflowingAdd
            + OverflowingSub
            + WrappingAdd
            + std::fmt::Display,
    > std::fmt::Display for MontgomeryModint<M, T>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.val())
    }
}

#[cfg(test)]
mod tests {
    use super::{Mod998244353, MontgomeryModint};

    #[test]
    fn montgomery_modint_test() {
        type Modint = MontgomeryModint<Mod998244353>;

        assert_eq!(Modint::zero().val(), 0);
        assert_eq!(Modint::one().val(), 1);
        assert_eq!(Modint::new(10).val(), 10);

        const A: i64 = 347384953;
        const B: i64 = 847362948;
        let a = Modint::new(A);
        let b = Modint::new(B);
        assert_eq!((a + b).val(), 196503548);
        assert_eq!((a - b).val(), 498266358);
        assert_eq!((a * b).val(), 17486571);
        assert_eq!(a.pow(B).val(), 860108694);
        assert_eq!((a / b).val(), 748159151);
    }

    #[test]
    fn montgomery_modint_u32_test() {
        type Modint = MontgomeryModint<Mod998244353<u32>, u32>;

        assert_eq!(Modint::zero().val(), 0);
        assert_eq!(Modint::one().val(), 1);
        assert_eq!(Modint::new(10).val(), 10);

        const A: u32 = 347384953;
        const B: u32 = 847362948;
        let a = Modint::new(A);
        let b = Modint::new(B);
        assert_eq!((a + b).val(), 196503548);
        assert_eq!((a - b).val(), 498266358);
        assert_eq!((a * b).val(), 17486571);
        assert_eq!(a.pow(B).val(), 860108694);
        assert_eq!((a / b).val(), 748159151);
    }
}
