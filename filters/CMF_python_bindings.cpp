#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "CMF_library.h"

namespace py = pybind11;

/*
 * Python binding for medfiltCirc2D.
 *
 * Parameters:
 *   y: 2D NumPy array (double) representing the image.
 *   R: Filter parameter R (must be odd integer)
 *   T: Filter parameter T (must be odd integer)
 *
 * Returns:
 *   A new 2D NumPy array with the filtered result.
 */
py::array_t<double> medfiltCirc2D_wrapper(py::array_t<double> y, int R, int T) {
    py::buffer_info buf_y = y.request();
    if (buf_y.ndim != 2)
        throw std::runtime_error("Input image y must be 2D");
    
    int M = static_cast<int>(buf_y.shape[0]);
    int N = static_cast<int>(buf_y.shape[1]);

    // Expected strides for a Fortran-contiguous array:
    size_t expected_stride0 = sizeof(double);
    size_t expected_stride1 = M * sizeof(double);
    if (buf_y.strides[0] != expected_stride0 || buf_y.strides[1] != expected_stride1) {
        throw std::runtime_error("Input array y must be Fortran contiguous (column-major). "
                                 "Please convert the array using np.asfortranarray() on the Python side.");
    }
    
    // Create output array with Fortran strides
    std::vector<size_t> shape = {static_cast<size_t>(M), static_cast<size_t>(N)};
    std::vector<size_t> fortran_strides = {sizeof(double), M * sizeof(double)};
    auto u = py::array_t<double>(shape, fortran_strides);
    
    // Get pointers
    py::buffer_info buf_u = u.request();
    double* ptr_y = static_cast<double*>(buf_y.ptr);
    double* ptr_u = static_cast<double*>(buf_u.ptr);
    
    // Call the underlying C++ function
    medfiltCirc2D(ptr_u, ptr_y, M, N, R, T);
    
    return u;
}

// Wrapper for the quantized circle median filter.
// y: 2D NumPy array (input image)
// R, T: filter parameters (must be odd integers)
// v: 1D NumPy array of quantization values.
py::array_t<double> medfiltCirc2DQuant_wrapper(py::array_t<double> y, int R, int T, py::array_t<double> v) {
    py::buffer_info buf_y = y.request();
    if (buf_y.ndim != 2)
        throw std::runtime_error("Input image y must be 2D");
    
    int M = static_cast<int>(buf_y.shape[0]);
    int N = static_cast<int>(buf_y.shape[1]);

    py::buffer_info buf_v = v.request();
    if (buf_v.ndim != 1)
        throw std::runtime_error("Input array v must be 1D");
    
    int S = static_cast<int>(buf_v.shape[0]);

    // Create output array with Fortran strides
    std::vector<size_t> shape = {static_cast<size_t>(M), static_cast<size_t>(N)};
    std::vector<size_t> fortran_strides = {sizeof(double), M * sizeof(double)};
    auto u = py::array_t<double>(shape, fortran_strides);
    py::buffer_info buf_u = u.request();

    double* ptr_y = static_cast<double*>(buf_y.ptr);
    double* ptr_v = static_cast<double*>(buf_v.ptr);
    double* ptr_u = static_cast<double*>(buf_u.ptr);

    // Call the underlying C++ function.
    medfiltCirc2DQuant(ptr_u, ptr_y, M, N, R, T, ptr_v, S);

    return u;
}

PYBIND11_MODULE(pycirclemedianfilter, m) {
    m.doc() = "Circle Median Filter module using C++ and pybind11";
    m.def("medfilt_circ2d", &medfiltCirc2D_wrapper,
          "Apply circle median filter on a 2D image (non-quantized data)",
          // py::arg("u"), py::arg("y"), py::arg("R"), py::arg("T"));
          py::arg("y"), py::arg("R"), py::arg("T"));
    m.def("medfilt_circ2d_quant", &medfiltCirc2DQuant_wrapper,
          "Apply circle median filter for quantized data on a 2D image",
          py::arg("y"), py::arg("R"), py::arg("T"), py::arg("v"));
}