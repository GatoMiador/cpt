/*
 * psim.h
 *
 *  Created on: 25 de jun. de 2023
 *      Author: marcos
 */

#ifndef CPT_H_
#define CPT_H_

#include <array>
#include <cmath>

/** Primitive for manipulating arrays. **/
template<class Z, std::size_t N> struct ArrayPrimitive : public std::array<Z, N> {
	/** Default contructor.
	 *
	 *     These arrays are always zeroed.
	**/
	ArrayPrimitive() = default;

	ArrayPrimitive(const Z & d) {
		for (std::size_t c=0; c<N; c++)
			(*this)[c] = d;
	};

	/** Copy contructor. **/
	ArrayPrimitive(const ArrayPrimitive & v) { *this = v; }

	/** Convert this array to another one. **/
	template<class X> explicit operator X() const noexcept {
		X b;
		for (std::size_t c=0; c<b.size(); c++)
			b[c] = (*this)[c];

		return b;
	}

	/** Array type. **/
	typedef Z value_type;

	typedef value_type * iterator;

	typedef const value_type * const_iterator;

	/** Returns an iterator of the start of array. **/
	inline constexpr iterator begin() noexcept { return iterator(std::array<Z, N>::begin() ); }

	/** Returns an iterator of the end of array. **/
	inline constexpr iterator end() noexcept { return iterator(std::array<Z, N>::end() ); }

	/** Returns a constant iterator of the start of the array. **/
	inline constexpr const_iterator begin() const noexcept { return const_iterator(std::array<Z, N>::begin() ); }

	/** Returns a constant iterator of the end of the array. **/
	inline constexpr const_iterator end() const noexcept { return const_iterator(std::array<Z, N>::end() ); }

	/** Compute a multiplication between two arrays. **/
	ArrayPrimitive operator*(const ArrayPrimitive & d) const noexcept {
		ArrayPrimitive r;
		for (std::size_t c=0; c<r.size(); c++)
			r[c] = (*this)[c] * d[c];

		return r;
	}

	/** Compute a division between two arrays. **/
	ArrayPrimitive operator/(const ArrayPrimitive & d) const noexcept {
		ArrayPrimitive r;
		for (std::size_t c=0; c<r.size(); c++)
			r[c] = (*this)[c] / d[c];

		return r;
	}

	/** Compute a division between an array and a numeric value. **/
	ArrayPrimitive operator/(const Z d) const noexcept {
		ArrayPrimitive r;
		for (std::size_t c=0; c<r.size(); c++)
			r[c] = (*this)[c] / d;

		return r;
	}

	/** Compute a subtraction between two arrays. **/
	ArrayPrimitive operator-(const ArrayPrimitive & d) const noexcept {
		ArrayPrimitive r;
		for (std::size_t c=0; c<r.size(); c++)
			r[c] = (*this)[c] - d[c];

	  return r;
	}

	/** Compute an addition between two arrays. **/
	ArrayPrimitive operator+(const ArrayPrimitive & d) const noexcept {
		ArrayPrimitive r;
		for (std::size_t c=0; c<r.size(); c++)
			r[c] = (*this)[c] + d[c];

	  return r;
	}

	/** Compute an addition between this array and another one. **/
	ArrayPrimitive& operator+=(const ArrayPrimitive& d) noexcept {
		for (std::size_t c=0; c<this->size(); c++)
			(*this)[c] += d[c];

		return *this;
	}

	/** Compute an subtraction between this array and another one. **/
	ArrayPrimitive& operator-=(const ArrayPrimitive& d) noexcept {
		for (std::size_t c=0; c<this->size(); c++)
			(*this)[c] -= d[c];

		return *this;
	}
};

template<unsigned int PHASES, unsigned int SAMPLE_RATE, unsigned int FREQ> class CPT {
public:
	/** Type for storing final values.
	 *
	 *     This array can be float or double depending on the tradeoff between
	 * performance and accuracy.
	**/
	typedef ArrayPrimitive<double, PHASES> power_vector;

	/** Type for storing values being computed. **/
	typedef ArrayPrimitive<double, PHASES> power_double_vector;

	/** Computes the square root of an array.
	 *
	 * @tparam T Type of the array
	 * @param v Array
	 * @return Result
	**/
	template<typename T> static ArrayPrimitive<T, PHASES> sqrtm(const ArrayPrimitive<T, PHASES> & v) noexcept {
		ArrayPrimitive<T, PHASES> r;
		for (std::size_t c=0; c<r.size(); c++)
			r[c] = sqrt(v[c]);

		return r;
	}

	/** Computes a Moving Average Filter (MAF). **/
	template<class C1 = power_vector, class C2 = power_double_vector>  class MAF {
	public:
		/** Inputs data and calculate the new MAF value
		 *
		 * @param v Input values
		 * @return This class for chaining purposes
		**/
		MAF & feed(const C1 & v) noexcept {
			for (std::size_t c=0; c<full.size(); c++) {
				full[c] -= data[index][c];
				full[c] += data[index][c] = v[c];
			}

			if (++index >= full.size() )
				index = 0;

			return *this;
		}

		/** Returns the last result of the MAF filter. **/
		C1 result(void) const noexcept {
			const auto a = full / data.size();

			return C1(a);
		}

		/** Returns the size of the MAF's internal array. **/
		constexpr std::size_t size() const noexcept {
			return data.size();
		};

	private:
		C2 full { 0 };
		std::array<C1, (SAMPLE_RATE/ FREQ)> data { 0 };
		std::size_t index = 0;
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
			const auto r = power_double_vector(v);

			val += r / SAMPLE_RATE;
			maf.feed(r);

			return *this;
		}

		/** Returns the last result of the MAF filter. **/
		power_vector result(void) const noexcept {
			const auto r = val - maf.result();

			return power_vector(r);
		}
	private:
		power_double_vector val { 0 };
		MAF<power_double_vector, power_double_vector> maf;
	};

	template<typename C> static void validate(C & v) {
		for (auto & o: v) {
			const auto c = std::fpclassify(o);
			if ( (c == FP_INFINITE) || (c == FP_NAN) )
				o = 0;
		}
	}

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
};

	struct Result feed(const power_vector & u, const power_vector & i) {
		struct Result r;

		r.u = u;
		r.i = i;

		// Compute the instantaeous active power per phase
		r.p = u * i;

		// Compute the average active power per phase
		r.P = p.feed(r.p).result();

		// Compute the unbiased integral of the voltage
		ui.feed(u);

		const auto _ui = ui.result();

		// Compute the instantaeous reactive energy per phase
		r.w = _ui * i;

		// Compute the average active power per phase
		r.W = w.feed(r.w).result();

		const auto sq_U = sq_u.feed(u * u).result();

		// Compute the RMS voltage
		r.U = sqrtm(sq_U);

		// Compute the instantaeous active current per phase
		r.ia = r.P * u / sq_U;

		// Validate ia because it is possible generation an invalid value
		validate(r.ia);

		const auto sq_Ui = sq_ui.feed(_ui * _ui).result();

		// Compute the instantaeous reactive current per phase
		r.ir = r.W * _ui / sq_Ui;

		// Validate ia because it is possible generation an invalid value
		validate(r.ir);

		// Compute instantaeous void current per phase.
		r.iv = i - r.ia - r.ir;

		return r;
	};
private:
	/** Stores the unbiased integral of the instantaneous voltage. **/
	UIntegral ui;

	/** Stores the moving average of the instantaneous active power. **/
	MAF<> p;

	/** Stores the moving average of the instantaneous active power. **/
	MAF<> w;

	/** Stores the square of the RMS of the voltage multiplied by the samples. **/
	MAF<> sq_u;

	/** Stores the square of the RMS of the unbiased integral of the voltage multiplied by the samples. **/
	MAF<> sq_ui;
};

#endif /* CPT_H_ */
