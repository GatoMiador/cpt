/*
 * psim.h
 *
 *  Created on: 25 de jun. de 2023
 *      Author: marcos
 */

#ifndef CPT_H_
#define CPT_H_

#include <math.h>

/** Calculates de CPT (Conservative Power Theory) parcels.
 *
 * https://ieeexplore.ieee.org/document/5638628
 *
 * @tparam PHASES Number of phases
 * @tparam SAMPLING_RATE Sampling rate
 * @tparam FREQ Mains frequency
**/
template<unsigned int PHASES, unsigned int SAMPLING_RATE, unsigned int FREQ> class CPT {
public:
	/** Class for manipulating arrays. **/
	template<typename Z, size_t N> struct Array {
		typedef Z value_type;

		typedef value_type * iterator;

		typedef const value_type * const_iterator;

		typedef size_t size_type;

		Array() = default;

		Array(const Z & d) {
			for (size_type c=0; c<N; c++)
				(*this)[c] = d;
		};

		/** Copy contructor. **/
		Array(const Array & v) { *this = v; }

		/** Convert this array to another one. **/
		template<typename X> explicit operator X() const noexcept {
			X b;
			for (size_type c=0; c<b.size(); c++)
				b[c] = (*this)[c];

			return b;
		}

		/** This is the [] operator. **/
		constexpr Z& operator[](size_type n) noexcept { return data[n]; }

		/** This is the const [] operator. **/
		constexpr const Z& operator[](size_type n) const noexcept { return data[n]; }

		/** Returns an iterator of the start of array. **/
		inline constexpr iterator begin() noexcept { return iterator(&data[0]); }

		/** Returns an iterator of the end of array. **/
		inline constexpr iterator end() noexcept { return iterator(&data[N]); }

		/** Returns a constant iterator of the start of the array. **/
		inline constexpr const_iterator begin() const noexcept { return const_iterator(&data[0]); }

		/** Returns a constant iterator of the end of the array. **/
		inline constexpr const_iterator end() const noexcept { return const_iterator(&data[N]); }

		constexpr size_type size() const noexcept { return N; }

		/** Computes a multiplication between two arrays. **/
		Array operator*(const Array & d) const noexcept {
			Array r;
			for (size_type c=0; c<r.size(); c++)
				r[c] = (*this)[c] * d[c];

			return r;
		}

		/** Computes a multiplication between an array and a numeric value. **/
		Array operator*(const Z d) const noexcept {
			Array r;
			for (size_type c=0; c<r.size(); c++)
				r[c] = (*this)[c] * d;

			return r;
		}

		/** Computes a division between two arrays. **/
		Array operator/(const Array & d) const noexcept {
			Array r;
			for (size_type c=0; c<r.size(); c++)
				r[c] = (*this)[c] / d[c];

			return r;
		}

		/** Computes a division between an array and a numeric value. **/
		Array operator/(const Z d) const noexcept {
			Array r;
			for (size_type c=0; c<r.size(); c++)
				r[c] = (*this)[c] / d;

			return r;
		}

		/** Computes a subtraction between two arrays. **/
		Array operator-(const Array & d) const noexcept {
			Array r;
			for (size_type c=0; c<r.size(); c++)
				r[c] = (*this)[c] - d[c];

		  return r;
		}

		/** Computes an addition between two arrays. **/
		Array operator+(const Array & d) const noexcept {
			Array r;
			for (size_type c=0; c<r.size(); c++)
				r[c] = (*this)[c] + d[c];

		  return r;
		}

		/** Computes an addition between this array and another one. **/
		Array& operator+=(const Array& d) noexcept {
			for (size_type c=0; c<this->size(); c++)
				(*this)[c] += d[c];

			return *this;
		}

		/** Computes an subtraction between this array and another one. **/
		Array& operator-=(const Array& d) noexcept {
			for (size_type c=0; c<this->size(); c++)
				(*this)[c] -= d[c];

			return *this;
		}

		/** Computes the summation of the values in this array.
		 *
		 * @tparam T Type of the resulting number.
		 * @return Result.
		**/
		template<typename T = Z> T sum(void) const noexcept {
			T res = 0;
			for (const auto & a: data)
				res += a;
			return res;
		}
	private:
		/** Stores the actual array data. **/
		Z data[N];
	};

	/** Type for output numbers.
	 *
	 *     This number can be float or double depending on the tradeoff between
	 * performance and accuracy.
	 * 	   A good choice may be int or long long for data acquired from ADC
	 * since float types do not behave well when they are added.
	**/
	typedef double output_number;

	/** Type for input numbers.
	 *
	 * 	   A good choice may be short int or int for data acquired from ADC
	 * since float types do not behave well when they are added.
	**/
	typedef double input_number;

	/** Type for storing final values. **/
	typedef Array<output_number, PHASES> power_vector;

	/** Type for storing values being computed. **/
	typedef Array<input_number, PHASES> power_double_vector;

	/** Computes the square root of an array.
	 *
	 * @tparam T Type of the array
	 * @param v Array
	 * @return Result
	**/
	template<typename T> static Array<T, PHASES> sqrtm(const Array<T, PHASES> & v) noexcept {
		Array<T, PHASES> r;
		for (size_t c=0; c<r.size(); c++)
			r[c] = sqrt(v[c]);

		return r;
	}

	/** Computes the Moving Average Filter (MAF).
	 *
	 *  FIXME - This class are known to cause numeric errors when using floating point types,
	 * alternative strategies are advised for these cases.
	 * Void power and currents should be zero for some tests and they are not.
	**/
	template<typename C1 = power_vector, typename C2 = power_double_vector>  class MAF {
	public:
		/** Inputs data and calculate the new MAF value
		 *
		 * @param v Input values
		 * @return This class for chaining purposes
		**/
		MAF & feed(const C1 & v) noexcept {
			for (size_t c=0; c<full.size(); c++) {
				full[c] -= data[index][c];
				full[c] += data[index][c] = v[c];
			}

			if (++index >= (sizeof(data)/sizeof(data[0]) ) )
				index = 0;

			return *this;
		}

		/** Returns the last result of the MAF filter. **/
		C1 result(void) const noexcept {
			const auto r = full / (sizeof(data)/sizeof(data[0]) );

			return C1(r);
		}
	private:
		C2 full { 0 };
		C1 data[SAMPLING_RATE/FREQ] { 0 };
		size_t index = 0;
	};

	/** Computes the unbiased integral. **/
	class UIntegral {
	public:
		/** Inputs data and calculate the new integral value
		 *
		 * @param v Input values
		 * @return This class for chaining purposes
		**/
		UIntegral & feed(const power_vector & v) noexcept {
			const auto r = decltype(val)(v);

			val += r + old_sample;
			old_sample = r;
			maf.feed(val);

			return *this;
		}

		/** Returns the last result of the MAF filter. **/
		power_vector result(void) const noexcept {
			const auto r = (val - maf.result() ) * 2 * M_PI * FREQ / (2 * SAMPLING_RATE);

			return power_vector(r);
		}
	private:
		power_double_vector old_sample { 0 };
		power_double_vector val { 0 };
		MAF<power_double_vector, power_double_vector> maf;
	};

	/** Validates a number and zeros invalid data.
	 *
	 * @tparam C Type of the number.
	 * @param v Number to be validatd.
	**/
	template<typename C> static void validate(C & o) noexcept {
		const auto c = fpclassify(o);
		if ( (c == FP_INFINITE) || (c == FP_NAN) )
			o = 0;
	}

	/** Validates an array and zeros invalid data.
	 *
	 * @tparam C Type of array.
	 * @param v Array to be validatd.
	**/
	template<typename C> static void validate(Array<C, PHASES> & v) noexcept {
		for (auto & o: v)
			validate(o);
	}
private:

	/** Stores the unbiased integral (û) of the instantaneous voltage. **/
	UIntegral ui;

	/** Stores the moving average of the instantaneous active power. **/
	MAF<> p;

	/** Stores the moving average of the instantaneous active power. **/
	MAF<> w;

	/** Stores the square of the RMS ||U|| of the voltage multiplied by the samples. **/
	MAF<> sq_u;

	/** Stores the square of the RMS ||Û||  of the unbiased integral of the voltage multiplied by the samples. **/
	MAF<> sq_ui;

	/** Stores the square of the RMS of the active current. **/
	MAF<> sq_ia;

	/** Stores the square of the RMS of the reactive current. **/
	MAF<> sq_ir;

	/** Stores the square of the RMS of the void current. **/
	MAF<> sq_iv;

	/** Stores the square of the RMS of the balanced active current. **/
	MAF<> sq_iba;

	/** Stores the square of the RMS of the balanced active current. **/
	MAF<> sq_ibr;

	/** Stores the square of the RMS of the unbalanced active current. **/
	MAF<> sq_iua;

	/** Stores the square of the RMS of the unbalanced reactive current. **/
	MAF<> sq_iur;
public:

	/** Results of CPT calculation. **/
	struct Result {
		/** Instantaeous active power per phase. **/
		power_vector p;

		/** Average active power per phase. **/
		power_vector P;

		/** Instantaeous reactive energy per phase. **/
		power_vector w;

		/** Average reactive energy per phase. **/
		power_vector W;

		/** Instantaeous voltage per phase. **/
		power_vector u;

		/** RMS of the voltage. **/
		power_vector U;

		/** Instantaeous current per phase. **/
		power_vector i;

		/** Instantaeous active current per phase. **/
		power_vector ia;

		/** Instantaeous reactive current per phase. **/
		power_vector ir;

		/** Instantaeous void current per phase. **/
		power_vector iv;

		/** Instantaeous balanced active current per phase. **/
		power_vector iba;

		/** Instantaeous balanced reactive current per phase. **/
		power_vector ibr;

		/** Instantaeous unbalanced active current per phase. **/
		power_vector iua;

		/** Instantaeous unbalanced reactive current per phase. **/
		power_vector iur;

		struct Totals {
			/** Overall active power. **/
			output_number P;

			/** Overall reactive power. **/
			output_number Q;

			/** Overall unbalance power. **/
			output_number N;

			/** Overall void power. **/
			output_number V;

			/** Overall apparent power. **/
			output_number A;

			/** Reactivity factor. **/
			output_number rf;

			/** Unbalance factor. **/
			output_number uf;

			/** Non linearity factor. **/
			output_number df;

			/** Power factor. **/
			output_number pf;
		} t;
	};

	struct Result feed(const power_vector & u, const power_vector & i) noexcept {
		struct Result r;

		// Just copy these values
		r.u = u;
		r.i = i;

		// Compute the instantaeous active power per phase
		r.p = u * i;

		// Compute the average active power per phase
		r.P = p.feed(r.p).result();

		//! Compute the unbiased integral of the voltage
		const auto _ui = ui.feed(u).result();

		// Compute the instantaeous reactive energy per phase
		r.w = _ui * i;

		// Compute the average active power per phase
		r.W = w.feed(r.w).result();

		//! Compute mean square of voltage
		const auto sq_U = sq_u.feed(u * u).result();

		// Compute the RMS voltage
		r.U = sqrtm(sq_U);

		// Compute the instantaeous active current per phase
		r.ia = r.P * u / sq_U;

		// Validate ia because it is possible generation an invalid value
		validate(r.ia);

		//! Compute the mean square of the unbiased integral of the voltage
		const auto sq_Ui = sq_ui.feed(_ui * _ui).result();

		// Compute the instantaeous reactive current per phase
		r.ir = r.W * _ui / sq_Ui;

		// Validate ia because it is possible generation an invalid value
		validate(r.ir);

		// Compute instantaeous void current per phase.
		/* FIXME - These subtractions are known to cause numeric errors when using floating point types,
		 * alternative strategies are advised for these cases.
		 * Void power and currents should be zero for some tests and they are not.
		*/
		r.iv = i - r.ia - r.ir;

		//! Compute the total mean square of voltage
		const auto sq_UT = sq_U.sum();

		//! Compute the total root mean square of voltage
		auto UT = sqrt(sq_UT);

		validate(UT);

		//! Compute the mean square of the void current
		const auto sq_iV = sq_iv.feed(r.iv * r.iv).result();

		// Compute overall void power.
		r.t.V = UT * sqrt(sq_iV.sum() );

		validate(r.t.V);

		if (PHASES > 1) {
			// Compute instantaeous balanced active current per phase
			r.iba = u * r.P.sum() / sq_UT;

			validate(r.iba);

			//! Compute the total mean square of the unbiased integral of voltage
			const auto sq_UiT = sq_Ui.sum();

			// Compute instantaeous balanced reactive current per phase
			r.ibr = _ui * r.W.sum() / sq_UiT;

			validate(r.ibr);

			// Compute instantaeous unbalanced active current per phase.
			r.iua = r.ia - r.iba;

			// Compute instantaeous unbalanced reactive current per phase.
			r.iur = r.ir - r.ibr;

			// Compute overall active power.
			r.t.P = UT * sqrt(sq_iba.feed(r.iba * r.iba).result().sum() );

			validate(r.t.P);

			// Compute overall reactive power.
			r.t.Q = UT * sqrt(sq_ibr.feed(r.ibr * r.ibr).result().sum() );

			validate(r.t.Q);

			const auto sq_iUa = sq_iua.feed(r.iua * r.iua).result();
			const auto sq_iUr = sq_iur.feed(r.iur * r.iur).result();

			// Compute overall unbalance power.
			r.t.N = UT * sqrt( (sq_iUa + sq_iUr).sum() );

			validate(r.t.N);

			// Compute overall apparent power.
			r.t.A = sqrt(r.t.P*r.t.P + r.t.Q*r.t.Q + r.t.N*r.t.N + r.t.V*r.t.V);

			validate(r.t.A);

			// Compute unbalance factor.
			r.t.uf = r.t.N / sqrt(r.t.P*r.t.P + r.t.Q*r.t.Q + r.t.N*r.t.N);

			validate(r.t.uf);
		} else {
			// Compute overall active power.
			r.t.P = UT * sqrt(sq_ia.feed(r.ia * r.ia).result().sum() );

			validate(r.t.P);

			// Compute overall reactive power.
			r.t.Q = UT * sqrt(sq_ir.feed(r.ir * r.ir).result().sum() );

			validate(r.t.Q);

			// Compute overall apparent power.
			r.t.A = sqrt(r.t.P*r.t.P + r.t.Q*r.t.Q + r.t.V*r.t.V);

			validate(r.t.A);
		}

		// Compute reactivity factor.
		r.t.rf = r.t.Q / sqrt(r.t.P*r.t.P + r.t.Q*r.t.Q);

		validate(r.t.rf);

		// Compute non linearity factor.
		r.t.df = r.t.V / r.t.A;

		validate(r.t.df);

		// Compute power factor.
		r.t.pf = r.t.P / r.t.A;

		validate(r.t.pf);

		return r;
	};
};

#endif /* CPT_H_ */
