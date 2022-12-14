use super::common::{
    bit_reverse, radix_4_inner_complex, radix_4_inner_montgomery_modint, radix_4_inv_inner_complex,
    radix_4_inv_inner_montgomery_modint, radix_8_inner_complex, radix_8_inner_montgomery_modint,
    radix_8_inv_inner_complex, radix_8_inv_inner_montgomery_modint,
};
use super::fft_cache::FftCache;
use super::modint::{Mod998244353, MontgomeryModint};
use num::Complex;
use std::ops::{Add, Mul, MulAssign, Sub};

#[inline]
fn radix_2_kernel<T>(deg: usize, width: usize, root: T, a: &mut [T])
where
    T: Clone + Copy + Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T> + MulAssign,
{
    let offset = width >> 1;
    for top in (0..deg).step_by(width) {
        let (c0, c1) = (a[top], a[top + offset]);
        let (c0, c1) = (c0 + c1, c0 - c1);
        a[top] = c0;
        a[top + offset] = c1;
        let mut w = root;
        for base in top + 1..top + offset {
            let (c0, c1) = (a[base], a[base + offset] * w);
            let (w0, w1) = (c0 + c1, c0 - c1);
            a[base] = w0;
            a[base + offset] = w1;
            w *= root;
        }
    }
}

type Radix4Inner<T> = fn(T, T, T, T, &FftCache<T>) -> (T, T, T, T);
type Radix8Inner<T> = fn(T, T, T, T, T, T, T, T, &FftCache<T>) -> (T, T, T, T, T, T, T, T);

#[inline]
fn radix_4_kernel<T>(
    deg: usize,
    width: usize,
    a: &mut [T],
    cache: &FftCache<T>,
    twiddle: &[T],
    radix4_inner: Radix4Inner<T>,
) where
    T: Clone + Copy + Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T> + MulAssign,
{
    let log = width.trailing_zeros();
    let exp = deg >> log;
    let exp = [0, exp, exp * 2, exp * 3];
    let offset = width >> 2;
    for top in (0..deg).step_by(width) {
        let (id0, id1, id2, id3) = (top, top + offset, top + offset * 2, top + offset * 3);
        let (c0, c1, c2, c3) = (a[id0], a[id2], a[id1], a[id3]);
        let (c0, c1, c2, c3) = radix4_inner(c0, c1, c2, c3, cache);
        a[id0] = c0;
        a[id1] = c1;
        a[id2] = c2;
        a[id3] = c3;
        for base in top + 1..top + offset {
            let w1 = twiddle[(exp[1] * (base - top)) & (deg - 1)];
            let w2 = twiddle[(exp[2] * (base - top)) & (deg - 1)];
            let w3 = twiddle[(exp[3] * (base - top)) & (deg - 1)];

            let (id0, id1, id2, id3) = (base, base + offset, base + offset * 2, base + offset * 3);
            let (c0, c1, c2, c3) = (a[id0], a[id2] * w1, a[id1] * w2, a[id3] * w3);

            let (c0, c1, c2, c3) = radix4_inner(c0, c1, c2, c3, cache);

            a[id0] = c0;
            a[id1] = c1;
            a[id2] = c2;
            a[id3] = c3;
        }
    }
}

#[inline]
fn radix_8_kernel<T>(
    deg: usize,
    width: usize,
    a: &mut [T],
    cache: &FftCache<T>,
    twiddle: &[T],
    radix_8_inner: Radix8Inner<T>,
) where
    T: Clone + Copy + Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T> + MulAssign,
{
    let log = width.trailing_zeros();
    let exp = deg >> log;
    let exp = [0, exp, exp * 2, exp * 3, exp * 4, exp * 5, exp * 6, exp * 7];
    let offset = width >> 3;
    for top in (0..deg).step_by(width) {
        let (c0, c1, c2, c3, c4, c5, c6, c7) = radix_8_inner(
            a[top],
            a[top + offset * 4],
            a[top + offset * 2],
            a[top + offset * 6],
            a[top + offset],
            a[top + offset * 5],
            a[top + offset * 3],
            a[top + offset * 7],
            cache,
        );
        a[top] = c0;
        a[top + offset] = c1;
        a[top + offset * 2] = c2;
        a[top + offset * 3] = c3;
        a[top + offset * 4] = c4;
        a[top + offset * 5] = c5;
        a[top + offset * 6] = c6;
        a[top + offset * 7] = c7;
        for base in top + 1..top + offset {
            let w1 = twiddle[(exp[1] * (base - top)) & (deg - 1)];
            let w2 = twiddle[(exp[2] * (base - top)) & (deg - 1)];
            let w3 = twiddle[(exp[3] * (base - top)) & (deg - 1)];
            let w4 = twiddle[(exp[4] * (base - top)) & (deg - 1)];
            let w5 = twiddle[(exp[5] * (base - top)) & (deg - 1)];
            let w6 = twiddle[(exp[6] * (base - top)) & (deg - 1)];
            let w7 = twiddle[(exp[7] * (base - top)) & (deg - 1)];

            let (id0, id1, id2, id3, id4, id5, id6, id7) = (
                base,
                base + offset,
                base + offset * 2,
                base + offset * 3,
                base + offset * 4,
                base + offset * 5,
                base + offset * 6,
                base + offset * 7,
            );
            let (c0, c1, c2, c3, c4, c5, c6, c7) = (
                a[id0],
                a[id4] * w1,
                a[id2] * w2,
                a[id6] * w3,
                a[id1] * w4,
                a[id5] * w5,
                a[id3] * w6,
                a[id7] * w7,
            );

            let (c0, c1, c2, c3, c4, c5, c6, c7) =
                radix_8_inner(c0, c1, c2, c3, c4, c5, c6, c7, cache);

            a[id0] = c0;
            a[id1] = c1;
            a[id2] = c2;
            a[id3] = c3;
            a[id4] = c4;
            a[id5] = c5;
            a[id6] = c6;
            a[id7] = c7;
        }
    }
}

macro_rules! impl_butterfly {
    ( $t:tt, $radix2:ident, $radix2_inv:ident, $radix4:ident, $radix4_inv:ident, $radix8:ident, $radix8_inv:ident, $radix4_inner:ident, $radix4_inner_inv:ident, $radix8_inner:ident, $radix8_inner_inv:ident ) => {
        #[inline]
        pub fn $radix2(deg: usize, log: usize, a: &mut [$t], cache: &FftCache<$t>) {
            for i in 0..log {
                let width = 1 << (i + 1);
                let root = cache.prim_root(i + 1);
                radix_2_kernel(deg, width, root, a);
            }
        }

        #[inline]
        pub fn $radix2_inv(deg: usize, log: usize, a: &mut [$t], cache: &FftCache<$t>) {
            for i in 0..log {
                let width = 1 << (i + 1);
                let root = cache.prim_root_inv(i + 1);
                radix_2_kernel(deg, width, root, a);
            }
        }

        #[inline]
        pub fn $radix4(deg: usize, log: usize, a: &mut [$t], cache: &FftCache<$t>) {
            for i in (0..log).step_by(2) {
                if i + 1 == log {
                    let width = 1 << (i + 1);
                    let root = cache.prim_root(i + 1);
                    radix_2_kernel(deg, width, root, a);
                } else {
                    let width = 1 << (i + 2);
                    radix_4_kernel(deg, width, a, cache, cache.twiddle_factors(), $radix4_inner);
                }
            }
        }

        #[inline]
        pub fn $radix4_inv(deg: usize, log: usize, a: &mut [$t], cache: &FftCache<$t>) {
            let twiddle_factors = cache.twiddle_factors_inv();
            for i in (0..log).step_by(2) {
                if i + 1 == log {
                    let width = 1 << (i + 1);
                    let root = cache.prim_root_inv(i + 1);
                    radix_2_kernel(deg, width, root, a);
                } else {
                    let width = 1 << (i + 2);
                    radix_4_kernel(deg, width, a, cache, twiddle_factors, $radix4_inner_inv);
                }
            }
        }

        #[inline]
        pub fn $radix8(deg: usize, log: usize, a: &mut [$t], cache: &FftCache<$t>) {
            for i in (0..log).step_by(3) {
                if i + 1 == log {
                    let width = 1 << (i + 1);
                    let root = cache.prim_root(i + 1);
                    radix_2_kernel(deg, width, root, a);
                } else if i + 2 == log {
                    let width = 1 << (i + 2);
                    radix_4_kernel(deg, width, a, cache, cache.twiddle_factors(), $radix4_inner);
                } else {
                    let width = 1 << (i + 3);
                    radix_8_kernel(deg, width, a, cache, cache.twiddle_factors(), $radix8_inner);
                }
            }
        }

        #[inline]
        pub fn $radix8_inv(deg: usize, log: usize, a: &mut [$t], cache: &FftCache<$t>) {
            let twiddle_factors = cache.twiddle_factors_inv();
            for i in (0..log).step_by(3) {
                if i + 1 == log {
                    let width = 1 << (i + 1);
                    let root = cache.prim_root_inv(i + 1);
                    radix_2_kernel(deg, width, root, a);
                } else if i + 2 == log {
                    let width = 1 << (i + 2);
                    radix_4_kernel(deg, width, a, cache, twiddle_factors, $radix4_inner_inv);
                } else {
                    let width = 1 << (i + 3);
                    radix_8_kernel(deg, width, a, cache, twiddle_factors, $radix8_inner_inv);
                }
            }
        }
    };
}

type Complexf32 = Complex<f32>;
type Complexf64 = Complex<f64>;
impl_butterfly!(
    Complexf32,
    cooley_tukey_radix_2_butterfly_f32,
    cooley_tukey_radix_2_butterfly_inv_f32,
    cooley_tukey_radix_4_butterfly_f32,
    cooley_tukey_radix_4_butterfly_inv_f32,
    cooley_tukey_radix_8_butterfly_f32,
    cooley_tukey_radix_8_butterfly_inv_f32,
    radix_4_inner_complex,
    radix_4_inv_inner_complex,
    radix_8_inner_complex,
    radix_8_inv_inner_complex
);
impl_butterfly!(
    Complexf64,
    cooley_tukey_radix_2_butterfly,
    cooley_tukey_radix_2_butterfly_inv,
    cooley_tukey_radix_4_butterfly,
    cooley_tukey_radix_4_butterfly_inv,
    cooley_tukey_radix_8_butterfly,
    cooley_tukey_radix_8_butterfly_inv,
    radix_4_inner_complex,
    radix_4_inv_inner_complex,
    radix_8_inner_complex,
    radix_8_inv_inner_complex
);

type MontMint998244353 = MontgomeryModint<Mod998244353<u32>, u32>;
impl_butterfly!(
    MontMint998244353,
    cooley_tukey_radix_2_butterfly_montgomery_modint,
    cooley_tukey_radix_2_butterfly_inv_montgomery_modint,
    cooley_tukey_radix_4_butterfly_montgomery_modint,
    cooley_tukey_radix_4_butterfly_inv_montgomery_modint,
    cooley_tukey_radix_8_butterfly_montgomery_modint,
    cooley_tukey_radix_8_butterfly_inv_montgomery_modint,
    radix_4_inner_montgomery_modint,
    radix_4_inv_inner_montgomery_modint,
    radix_8_inner_montgomery_modint,
    radix_8_inv_inner_montgomery_modint
);

macro_rules! impl_fft_inner {
    ( $t:ty, $butterfly:ident, $deg:ident, $log:ident, $a:ident, $cache:ident, $epilogue:expr ) => {{
        $butterfly($deg, $log, $a, &($cache));
        $epilogue($deg, $a);
    }};
}

macro_rules! impl_fft_template {
    ( $t:ty, $fname:ident, $butterfly:ident, $inner:expr ) => {
        pub fn $fname(a: &mut [$t]) {
            let deg = a.len();
            let log = deg.trailing_zeros() as usize;
            debug_assert_eq!(a.len(), 1 << log);
            bit_reverse(deg, log, a);
            let cache = FftCache::<$t>::new(log);
            impl_fft_inner!($t, $butterfly, deg, log, a, cache, $inner)
        }
    };
}

macro_rules! impl_fft {
    ( $t:ty, $fname:ident, $butterfly:ident, $fname_inv:ident, $butterfly_inv:ident, $inner_inv:expr) => {
        impl_fft_template!($t, $fname, $butterfly, {});
        impl_fft_template!($t, $fname_inv, $butterfly_inv, $inner_inv);
    };
}

macro_rules! impl_fft_all_radix {
    ( $t:ty, $rad2:ident, $butterfly2:ident, $rad2inv:ident, $butterfly2inv:ident, $rad4:ident, $butterfly4:ident, $rad4inv:ident, $butterfly4inv:ident, $rad8:ident, $butterfly8:ident, $rad8inv:ident, $butterfly8inv:ident, $inner_inv:expr) => {
        impl_fft!($t, $rad2, $butterfly2, $rad2inv, $butterfly2inv, $inner_inv);
        impl_fft!($t, $rad4, $butterfly4, $rad4inv, $butterfly4inv, $inner_inv);
        impl_fft!($t, $rad8, $butterfly8, $rad8inv, $butterfly8inv, $inner_inv);
    };
}

impl_fft_all_radix!(
    Complexf64,
    fft_cooley_tukey_radix_2,
    cooley_tukey_radix_2_butterfly,
    ifft_cooley_tukey_radix_2,
    cooley_tukey_radix_2_butterfly_inv,
    fft_cooley_tukey_radix_4,
    cooley_tukey_radix_4_butterfly,
    ifft_cooley_tukey_radix_4,
    cooley_tukey_radix_4_butterfly_inv,
    fft_cooley_tukey_radix_8,
    cooley_tukey_radix_8_butterfly,
    ifft_cooley_tukey_radix_8,
    cooley_tukey_radix_8_butterfly_inv,
    |deg: usize, a: &mut [Complexf64]| a.iter_mut().for_each(|c| *c = c.conj() / deg as f64)
);

impl_fft_all_radix!(
    Complexf32,
    fft_cooley_tukey_radix_2_f32,
    cooley_tukey_radix_2_butterfly_f32,
    ifft_cooley_tukey_radix_2_f32,
    cooley_tukey_radix_2_butterfly_inv_f32,
    fft_cooley_tukey_radix_4_f32,
    cooley_tukey_radix_4_butterfly_f32,
    ifft_cooley_tukey_radix_4_f32,
    cooley_tukey_radix_4_butterfly_inv_f32,
    fft_cooley_tukey_radix_8_f32,
    cooley_tukey_radix_8_butterfly_f32,
    ifft_cooley_tukey_radix_8_f32,
    cooley_tukey_radix_8_butterfly_inv_f32,
    |deg: usize, a: &mut [Complexf32]| a.iter_mut().for_each(|c| *c = c.conj() / deg as f32)
);

impl_fft_all_radix!(
    MontMint998244353,
    fft_cooley_tukey_radix_2_montgomery_modint,
    cooley_tukey_radix_2_butterfly_montgomery_modint,
    ifft_cooley_tukey_radix_2_montgomery_modint,
    cooley_tukey_radix_2_butterfly_inv_montgomery_modint,
    fft_cooley_tukey_radix_4_montgomery_modint,
    cooley_tukey_radix_4_butterfly_montgomery_modint,
    ifft_cooley_tukey_radix_4_montgomery_modint,
    cooley_tukey_radix_4_butterfly_inv_montgomery_modint,
    fft_cooley_tukey_radix_8_montgomery_modint,
    cooley_tukey_radix_8_butterfly_montgomery_modint,
    ifft_cooley_tukey_radix_8_montgomery_modint,
    cooley_tukey_radix_8_butterfly_inv_montgomery_modint,
    |deg: usize, a: &mut [MontMint998244353]| {
        let inv = MontMint998244353::new(deg as u32).inv();
        a.iter_mut().for_each(|c| *c *= inv)
    }
);

#[cfg(test)]
mod tests {
    use super::super::modint::{Mod998244353, MontgomeryModint};
    use super::{
        fft_cooley_tukey_radix_2, fft_cooley_tukey_radix_4,
        fft_cooley_tukey_radix_4_montgomery_modint, fft_cooley_tukey_radix_8,
        fft_cooley_tukey_radix_8_montgomery_modint, ifft_cooley_tukey_radix_2,
        ifft_cooley_tukey_radix_4, ifft_cooley_tukey_radix_8,
        ifft_cooley_tukey_radix_8_montgomery_modint,
    };
    use num::Complex;

    fn calc_diff(a: &[Complex<f64>], b: &[Complex<f64>]) -> f64 {
        let mut diff_max = 0.0;
        for (d, v) in a.into_iter().zip(b.into_iter()) {
            if (d.re - v.re).abs() > diff_max {
                diff_max = (d.re - v.re).abs();
            }
            if (d.im - v.im).abs() > diff_max {
                diff_max = (d.im - v.im).abs();
            }
        }

        diff_max
    }

    #[test]
    fn cooley_tukey_radix_2_test() {
        let data: Vec<Complex<f64>> = (1..=16).map(|v| (v as f64).into()).collect();
        let mut data1 = data.clone();
        fft_cooley_tukey_radix_2(&mut data1);
        ifft_cooley_tukey_radix_2(&mut data1);
        let diff_max = calc_diff(&data, &data1);
        assert!(diff_max < 1e-10);
    }

    #[test]
    fn cooley_tukey_radix_4_test() {
        let data: Vec<Complex<f64>> = (1..=16).map(|v| (v as f64).into()).collect();
        let mut data1 = data.clone();
        fft_cooley_tukey_radix_4(&mut data1);
        ifft_cooley_tukey_radix_4(&mut data1);
        let diff_max = calc_diff(&data, &data1);
        assert!(diff_max < 1e-10);
    }

    #[test]
    fn cooley_tukey_radix_8_test() {
        let data: Vec<Complex<f64>> = (1..=16).map(|v| (v as f64).into()).collect();
        let mut data1 = data.clone();
        fft_cooley_tukey_radix_8(&mut data1);
        ifft_cooley_tukey_radix_8(&mut data1);
        let diff_max = calc_diff(&data, &data1);
        assert!(diff_max < 1e-10);
    }

    #[test]
    fn cooley_tukey_radix_8_montgomery_modint_test() {
        type MontMint998244353 = MontgomeryModint<Mod998244353<u32>, u32>;
        let data: Vec<MontMint998244353> = (1..=16).map(|v| MontMint998244353::new(v)).collect();
        let mut data1 = data.clone();
        fft_cooley_tukey_radix_8_montgomery_modint(&mut data1);
        ifft_cooley_tukey_radix_8_montgomery_modint(&mut data1);
        assert_eq!(data, data1);
    }

    #[test]
    fn cooley_tukey_radix_2_radix_4_compare_test() {
        let data: Vec<Complex<f64>> = (1..=16).map(|v| (v as f64).into()).collect();
        let mut data1 = data.clone();
        let mut data2 = data.clone();
        fft_cooley_tukey_radix_2(&mut data1);
        fft_cooley_tukey_radix_4(&mut data2);
        let diff_max = calc_diff(&data1, &data2);
        assert!(diff_max < 1e-10);

        ifft_cooley_tukey_radix_2(&mut data1);
        let diff_max = calc_diff(&data, &data1);
        assert!(diff_max < 1e-10);

        ifft_cooley_tukey_radix_4(&mut data2);
        let diff_max = calc_diff(&data, &data2);
        assert!(diff_max < 1e-10);
    }

    #[test]
    fn gentleman_sande_radix_4_radix_8_compare_test() {
        type MontMint998244353 = MontgomeryModint<Mod998244353<u32>, u32>;
        let data: Vec<MontMint998244353> = (1..=16).map(|v| MontMint998244353::new(v)).collect();
        let mut data1 = data.clone();
        let mut data2 = data.clone();
        fft_cooley_tukey_radix_8_montgomery_modint(&mut data1);
        fft_cooley_tukey_radix_4_montgomery_modint(&mut data2);
        assert_eq!(data1, data2);
    }
}
