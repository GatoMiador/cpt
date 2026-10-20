/*
 * _cpt.cpp
 *
 *  Python bindings (nanobind) for the C++ CPT implementation in cpt.h.
 *
 *  The CPT class is a template resolved at compile time, so a fixed set of
 *  configurations (phases x sampling rate x mains frequency x number type)
 *  is instantiated here. To add a configuration, add a line to the list in
 *  NB_MODULE at the end of this file.
 */

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "../cpt.h"

namespace nb = nanobind;
using namespace nb::literals;

/** Result for single phase systems.
 *
 *     It leaves out the quantities defined only for polyphase systems
 * (iba, ibr, iua, iur, N and uf). CPT::feed() does not compute them when
 * PHASES == 1, so they would hold uninitialized values with the default
 * CPTResult.
**/
template<typename power_vector, typename output_vector>
struct MonoResult {
	power_vector p, P, w, W, u, U, i, ia, ir, iv;

	struct Totals {
		output_vector P, Q, V, A, rf, df, pf;
	} t;
};

/** Selects the result type for the number of phases. **/
template<unsigned int PHASES> struct ResultFor {
	template<typename A, typename B> using type = CPTconcepts::CPTResult<A, B>;
};

template<> struct ResultFor<1> {
	template<typename A, typename B> using type = MonoResult<A, B>;
};

/** Names of the per phase quantities, in output order. **/
static constexpr std::string_view phase_names[] = {
	"u", "i", "p", "P", "w", "W", "U", "ia", "ir", "iv", "iba", "ibr", "iua", "iur"
};

/** Names of the overall quantities, in output order (prefixed with "t." in 'fields'). **/
static constexpr std::string_view total_names[] = {
	"P", "Q", "N", "V", "A", "rf", "uf", "df", "pf"
};

constexpr std::size_t n_phase = std::size(phase_names);
constexpr std::size_t n_total = std::size(total_names);

/** Returns a pointer to the per phase quantity 'k' of a result, or nullptr if the result does not have it. **/
template<typename R> static const auto * phase_ptr(const R & r, std::size_t k) noexcept {
	using V = decltype(r.u);
	const V * p = nullptr;
	switch (k) {
	case 0:  p = &r.u; break;
	case 1:  p = &r.i; break;
	case 2:  p = &r.p; break;
	case 3:  p = &r.P; break;
	case 4:  p = &r.w; break;
	case 5:  p = &r.W; break;
	case 6:  p = &r.U; break;
	case 7:  p = &r.ia; break;
	case 8:  p = &r.ir; break;
	case 9:  p = &r.iv; break;
	case 10: if constexpr (requires { r.iba; }) p = &r.iba; break;
	case 11: if constexpr (requires { r.ibr; }) p = &r.ibr; break;
	case 12: if constexpr (requires { r.iua; }) p = &r.iua; break;
	case 13: if constexpr (requires { r.iur; }) p = &r.iur; break;
	}
	return p;
}

/** Returns a pointer to the overall quantity 'k' of a result, or nullptr if the result does not have it. **/
template<typename R> static const auto * total_ptr(const R & r, std::size_t k) noexcept {
	using V = decltype(r.t.P);
	const V * p = nullptr;
	switch (k) {
	case 0: p = &r.t.P; break;
	case 1: p = &r.t.Q; break;
	case 2: if constexpr (requires { r.t.N; }) p = &r.t.N; break;
	case 3: p = &r.t.V; break;
	case 4: p = &r.t.A; break;
	case 5: p = &r.t.rf; break;
	case 6: if constexpr (requires { r.t.uf; }) p = &r.t.uf; break;
	case 7: p = &r.t.df; break;
	case 8: p = &r.t.pf; break;
	}
	return p;
}

/** Creates a NumPy array that owns 'data' (allocated with new[]). **/
template<typename X> static nb::object to_numpy(X * data, std::initializer_list<std::size_t> shape) {
	nb::capsule owner(data, [](void * p) noexcept { delete[] static_cast<X *>(p); });
	return nb::ndarray<nb::numpy, X>(data, shape, owner).cast();
}

/** Python facing wrapper of one CPT configuration. **/
template<unsigned int PHASES, unsigned int SAMPLING_RATE, unsigned int FREQ, typename T>
class Wrapper {
public:
	using Core = CPT<PHASES, SAMPLING_RATE, FREQ, 1, T, T, ResultFor<PHASES>::template type>;

	using Input = nb::ndarray<const T, nb::ndim<2>, nb::c_contig, nb::device::cpu>;

	Wrapper() : core(std::make_unique<Core>()) { }

	/** Clears the internal state (moving average windows and integrator). **/
	void reset() {
		core = std::make_unique<Core>();
		count = 0;
	}

	/** Number of samples fed since construction or the last reset. **/
	std::uint64_t samples() const noexcept { return count; }

	/** Feeds a block of samples.
	 *
	 * @param u Voltages, shape (N, PHASES)
	 * @param i Currents, shape (N, PHASES)
	 * @param fields Quantities to return; None returns all the configuration computes
	 * @param step Keeps one result every 'step' samples (counted from the first sample fed)
	 * @return Dictionary with the selected quantities and 'n', the sample number of each row
	**/
	nb::dict feed_block(Input u, Input i, std::optional<std::vector<std::string>> fields, std::size_t step) {
		if ( (u.shape(1) != PHASES) || (i.shape(1) != PHASES) )
			throw std::invalid_argument("u and i must have shape (N, " + std::to_string(PHASES) + ")");
		if (u.shape(0) != i.shape(0) )
			throw std::invalid_argument("u and i must have the same number of samples");
		if (step == 0)
			throw std::invalid_argument("step must be at least 1");

		// Selects the quantities
		bool want_phase[n_phase] = { };
		bool want_total[n_total] = { };
		select(fields, want_phase, want_total);

		const std::size_t n = u.shape(0);
		const std::uint64_t first = count;
		const std::size_t rows = (first + n) / step - first / step;

		// Allocates the outputs before releasing the GIL
		T * phase_data[n_phase] = { };
		T * total_data[n_total] = { };
		auto * index = new std::int64_t[rows];
		for (std::size_t k = 0; k < n_phase; k++)
			if (want_phase[k])
				phase_data[k] = new T[rows * PHASES];
		for (std::size_t k = 0; k < n_total; k++)
			if (want_total[k])
				total_data[k] = new T[rows];

		{
			nb::gil_scoped_release release;

			const T * pu = u.data();
			const T * pi = i.data();
			std::size_t row = 0;
			typename Core::power_vector vu, vi;

			for (std::size_t s = 0; s < n; s++) {
				for (unsigned int c = 0; c < PHASES; c++) {
					vu[c] = pu[s * PHASES + c];
					vi[c] = pi[s * PHASES + c];
				}

				const auto r = core->feed(vu, vi);

				if (++count % step)
					continue;

				index[row] = static_cast<std::int64_t>(count - 1);
				for (std::size_t k = 0; k < n_phase; k++)
					if (phase_data[k]) {
						const auto * v = phase_ptr(r, k);
						for (unsigned int c = 0; c < PHASES; c++)
							phase_data[k][row * PHASES + c] = (*v)[c];
					}
				for (std::size_t k = 0; k < n_total; k++)
					if (total_data[k])
						total_data[k][row] = *total_ptr(r, k);
				row++;
			}
		}

		nb::dict out, totals;
		out["n"] = to_numpy(index, { rows });
		for (std::size_t k = 0; k < n_phase; k++)
			if (phase_data[k])
				out[nb::str(phase_names[k].data(), phase_names[k].size())] = to_numpy(phase_data[k], { rows, PHASES });
		for (std::size_t k = 0; k < n_total; k++)
			if (total_data[k])
				totals[nb::str(total_names[k].data(), total_names[k].size())] = to_numpy(total_data[k], { rows });
		if (totals.size() )
			out["totals"] = totals;

		return out;
	}

	/** Names of the quantities this configuration computes ("t." prefix for the overall ones). **/
	static std::vector<std::string> available() {
		std::vector<std::string> names;
		const typename Core::Result r { };
		for (std::size_t k = 0; k < n_phase; k++)
			if (phase_ptr(r, k) )
				names.emplace_back(phase_names[k]);
		for (std::size_t k = 0; k < n_total; k++)
			if (total_ptr(r, k) )
				names.push_back("t." + std::string(total_names[k]) );
		return names;
	}

private:
	/** Translates the field names into flags, validating them. **/
	static void select(const std::optional<std::vector<std::string>> & fields, bool * want_phase, bool * want_total) {
		const typename Core::Result r { };

		if (!fields) {
			for (std::size_t k = 0; k < n_phase; k++)
				want_phase[k] = phase_ptr(r, k) != nullptr;
			for (std::size_t k = 0; k < n_total; k++)
				want_total[k] = total_ptr(r, k) != nullptr;
			return;
		}

		for (const auto & f: *fields) {
			const bool is_total = f.starts_with("t.");
			const std::string_view name = is_total ? std::string_view(f).substr(2) : std::string_view(f);
			const auto * names = is_total ? total_names : phase_names;
			const std::size_t size = is_total ? n_total : n_phase;

			std::size_t k = 0;
			while ( (k < size) && (names[k] != name) )
				k++;

			if (k == size)
				throw std::invalid_argument("unknown field '" + f + "'");
			if ( (is_total ? (const void *) total_ptr(r, k) : (const void *) phase_ptr(r, k) ) == nullptr)
				throw std::invalid_argument("field '" + f + "' is not defined for single phase systems");

			(is_total ? want_total : want_phase)[k] = true;
		}
	}

	std::unique_ptr<Core> core;
	std::uint64_t count = 0;
};

/** Binds one configuration as class CPT_<phases>_<rate>_<freq>_<d|f>. **/
template<unsigned int PHASES, unsigned int SAMPLING_RATE, unsigned int FREQ, typename T>
static void bind(nb::module_ & m, nb::list & configurations) {
	static_assert(SAMPLING_RATE % FREQ == 0, "The sampling rate must be a multiple of the mains frequency");

	using W = Wrapper<PHASES, SAMPLING_RATE, FREQ, T>;
	const char * type = std::is_same_v<T, float> ? "float" : "double";
	const std::string name = "CPT_" + std::to_string(PHASES) + "_" + std::to_string(SAMPLING_RATE) + "_" +
		std::to_string(FREQ) + "_" + type[0];

	nb::class_<W>(m, name.c_str() )
		.def(nb::init<>() )
		.def("reset", &W::reset)
		.def_prop_ro("samples", &W::samples)
		.def("feed_block", &W::feed_block, "u"_a, "i"_a, "fields"_a = nb::none(), "step"_a = 1)
		.def_static("available", &W::available);

	configurations.append(nb::make_tuple(PHASES, SAMPLING_RATE, FREQ, type) );
}

/** Binds single and three phase, double and float, for one sampling rate and frequency. **/
template<unsigned int SAMPLING_RATE, unsigned int FREQ>
static void bind_rate(nb::module_ & m, nb::list & configurations) {
	bind<1, SAMPLING_RATE, FREQ, double>(m, configurations);
	bind<3, SAMPLING_RATE, FREQ, double>(m, configurations);
	bind<1, SAMPLING_RATE, FREQ, float>(m, configurations);
	bind<3, SAMPLING_RATE, FREQ, float>(m, configurations);
}

NB_MODULE(_cpt, m) {
	m.doc() = "C++ implementation of the Conservative Power Theory (cpt.h)";

	nb::list configurations;

	// 60 Hz
	bind_rate< 7680, 60>(m, configurations);	// 128 samples per cycle
	bind_rate<12000, 60>(m, configurations);	// 200 samples per cycle
	bind_rate<15360, 60>(m, configurations);	// 256 samples per cycle
	bind_rate<30720, 60>(m, configurations);	// 512 samples per cycle
	bind_rate<61440, 60>(m, configurations);	// 1024 samples per cycle

	// 50 Hz
	bind_rate< 6400, 50>(m, configurations);	// 128 samples per cycle
	bind_rate<12000, 50>(m, configurations);	// 240 samples per cycle
	bind_rate<12800, 50>(m, configurations);	// 256 samples per cycle
	bind_rate<25600, 50>(m, configurations);	// 512 samples per cycle
	bind_rate<51200, 50>(m, configurations);	// 1024 samples per cycle

	m.attr("configurations") = configurations;
}
