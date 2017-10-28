#!/usr/bin/env/ python3
"""
Objects for implementing spectral methods.

This module provides several methods contained within the SpectralToolbox class
for performing standard vector calculus operations in spectral space. It is
intended for use on a two-dimensional domain of equally-spaced gridpoints.

Classes
-------
- `SpectralToolbox` -- Superclass to implement standard techniques of spectral methods

Miscellaneous Functions for Unit Testing
----------------------------------------
- `test_func` -- Test function for unit testing
- `check_derivatives`
- `check_jacobian`
- `check_multiply`
- `check_laplacian`

All of the above compare the computed solution to an analytical solution.

Notes
-----
1) Invoking this module from the command line without any arguments will perform unit tests on the methods contained within.

2) The Fourier representation used is:

.. math:: u(x, y, t) = \\sum_{k_{1} = -N/2}^{N/2 -1} \\sum_{k_{2} = -N/2}^{N/2 - 1} \hat{u}(t) \\exp((2i\\pi/L)(k_{1}x + k_{2}y))

The spectrum is stored in the logical-to-mathematicians order, i.e. the wavenumbers are stored as:

+--------------------+--------------------+--------------------+--------------------+
|    (-N/2)(-N/2)    |  (-N/2)(-N/2 + 1)  |         ...        |   (-N/2)(N/2 - 1)  |
+--------------------+--------------------+--------------------+--------------------+
|  (-N/2 + 1)(-N/2)  |(-N/2 + 1)(-N/2 + 1)|         ...        |  (-N/2 + 1)(N/2 -1)|
+--------------------+--------------------+--------------------+--------------------+
|         ...        |         ...        |         ...        |         ...        |
+--------------------+--------------------+--------------------+--------------------+
|  (N/2 - 1)(-N/2)   | (N/2 - 1)(-N/2 + 1)|         ...        |  (N/2 - 1)(N/2 -1) |
+--------------------+--------------------+--------------------+--------------------+

With the Nyquist frequencies (corresponding to the most negative frequencies) located where one or both index = -N/2,
and the constant mode corresponding to :math:`k_{1} = 0`, :math:`k_{2} = 0` in the (N, N) position in the array.

2) The sign matrix, `sign_mat`, is used to relate the truncated Fourier series of a function, with wavenumbers that
run from [-N/2, N/2), with the FFT, with wavenumbers in [0, N). The point of doing this is to make the
operations in spectral space more intuitive from a mathematical perspective. We proceed by discretising the
desired function on an equispaced grid as above, noting that the sum runs from k = -N/2 to N/2 -1. In order to
make the sum correspond with the FFT, we perform the change of variables :math:`p = k + N/2`:

.. math:: f(x_{j}) = \\sum_{p=0}^{N-1} \\hat{f}(p-N/2) e^{2i\\pi k p/N}

.. math:: f(x_{j}) = \\sum_{p=0}^{N-1} \\hat{f}(p-N/2) (-1)^{p} e^{2i\\pi k/N}

Note that a factor of :math:`(-1)^{p}` is obtained. These alternating signs are implemented in the sign matrix
to permit a single multiplication operation.

See also
--------
numpy.fft, numpy

| Author: Adam G. Peddle
| Contact: ap553@exeter.ac.uk
| Version: 1.0
"""

import numpy as np

class SpectralToolbox:
    """
    Implements the required vector calculus techniques for spectral methods.

    This class stores the required constants as attributes to be used by its
    methods for spectral methods in 2 spatial dimensions.

    **Attributes**

    - `N` : The number of grid points in the domain, such that there are N X N equally-spaced points
    - `L` : The side length of the domain
    - `factor` : Factor arising in spectral differentiation on non-2pi domain (see below)
    - `sign_mat` : Matrix to reorder spectral space (see below)
    - `deriv_mat_x1`, deriv_mat_x2` : matrices for spectral differentiation in x1 and x2 directions (see calc_derivative method)
    - `laplace_op` : Array to implement Laplacian (see laplacian method)

    **Methods**

    - `dealias_pad` -- Pads the input array with zeros to prevent aliasing
    - `dealias_unpad` -- Removes the padding from dealias_pad
    - `forward_fft` -- Wrapper to perform Fast Fourier Transform with desired normalisation
    - `inverse_fft` -- Wrapper to invert FFT with desired normalisation
    - `calc_derivative` -- Computes 1st or 2nd order derivatives along x and/or y directions
    - `multiply_nonlinear` -- Multiplies two functions in a pseudo-spectral fashion
    - `solve_inverse_laplacian` -- Solves the inverse Laplacian problem with a chosen constant of integration
    - `jacobian1` -- Computes Jacobian of a single function
    - `jacobian` -- Computes Jacobian of two functions
    - `laplacian` -- Computes the Laplacian of a given function

    **Example**

    | ``>> st = SpectralToolbox(128, 2*np.pi)``
    | ``>> du_dx = st.calc_derivative(A, 'x')  # Compute x-derivative``
    """

    def __init__(self, N, L):
        self.N = N
        if N%2: # Methods are not currently implemented for odd values of N
            raise ValueError("N must be even")
        self.L = L
        self.factor = 2.0*np.pi/L

        # Create sign matrix
        self.sign_mat = np.zeros((2*self.N,2*self.N), dtype='complex')
        for k1 in range(0,2*self.N):
            for k2 in range(0,2*self.N):
                self.sign_mat[k1,k2] = (-1)**(k1+k2)

        # Create derivative matrices
        self.deriv_mat_x1 = np.zeros((self.N, self.N), dtype = complex)
        self.deriv_mat_x2 = np.zeros((self.N, self.N), dtype = complex)
        for k1 in range(-self.N//2,self.N//2):
            for k2 in range(-self.N//2,self.N//2):
                self.deriv_mat_x1[k1+self.N//2,k2+self.N//2] = 1j * k1 * self.factor
                self.deriv_mat_x2[k1+self.N//2,k2+self.N//2] = 1j * k2 * self.factor

        # Create laplacian matrix. Factor of multiplication is included
        # implicitly due to inclusion in derivative matrices.
        self.laplace_op = np.zeros((self.N, self.N), dtype = complex)
        self.laplace_op = self.deriv_mat_x1**2 + self.deriv_mat_x2**2

    def dealias_pad(self, A):
        """
        Pads the spectrum with zeros in high frequencies for aliasing control.

        This method doubles the size of the spectrum and sets all high frequencies to zero
        for aliasing control. This causes the resolution in real space (after inverse FFT)
        to be doubled. The higher frequencies are then discarded upon returning to
        spectral space. The array is doubled in a symmetric fashion.

        **Parameters**

        - `A` : the spectrum to be padded

        **Returns**

        - `B` : the padded spectrum

        **See also**

        dealias_unpad

        **Example**

        | ``>> A = np.ones((2,2), dtype = complex)``
        | ``>> st.dealias_pad(A)``
        | ``>> array([[ 0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j],``
        |          ``[ 0.+0.j,  1.+0.j,  1.+0.j,  0.+0.j],``
        |          ``[ 0.+0.j,  1.+0.j,  1.+0.j,  0.+0.j],``
        |          ``[ 0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j]])``

        """
        m, n = A.shape
        B = np.zeros((2*m, 2*n), dtype = complex)
        B[m//2:3*m//2, n//2:3*n//2] = A[:,:]
        return B

    def dealias_unpad(self, A):
        """
        Unpads the spectrum, removing high frequencies, for aliasing control.

        This method reverses the padding set up in the dealias_pad method,
        returning the size of the spectrum to the initial value upon returning
        to spectral space. The highest half of the frequencies are discarded.

        **Parameters**

        - `A` : the spectrum to be unpadded

        **Returns**

        - `B` : the unpadded spectrum

        **See also**

        dealias_pad

        **Example**

        | ``>> A = array([[ 0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j],``
        |            ``[ 0.+0.j,  1.+0.j,  1.+0.j,  0.+0.j],``
        |            ``[ 0.+0.j,  1.+0.j,  1.+0.j,  0.+0.j],``
        |            ``[ 0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j]])``
        | ``>> st.dealias_unpad(A)``
        | ``>> array([[ 1.+0.j,  1.+0.j],``
        |        ``[ 1.+0.j,  1.+0.j]])``
        """
        m, n = A.shape
        B = np.zeros((m//2, n//2), dtype = complex)
        B = A[m//4:3*m//4, n//4:3*n//4]
        return B

    def forward_fft(self, array_in):
        """
        Wrapper to implement the forward Fast Fourier Transform.

        This method implements the forward FFT (i.e. from real-space
        to spectral space). The values are normalised such that the spectral
        coefficient corresponds to a wave of exactly that amplitude, rather
        than the normalised value returned by the standard FFT. The sign matrix
        is used to yield a spectrum running from -N/2 to N/2 -1.

        **Parameters**

        - `array_in` : the array of real values

        **Returns**

        - `out` : the spectral coefficients, in the chosen framework

        **See Also**

        Note 2 in the class header (`sign matrix`), inverse_fft

        """
        # Find side length, as real array may or may not be doubled for
        # aliasing control
        side = array_in.shape[0]
        div_side = 1.0/side**2

        out = np.fft.fft2(self.sign_mat[0:side, 0:side]*array_in)*div_side
        return out

    def inverse_fft(self, array_in):
        """
        Wrapper to implement the inverse Fast Fourier Transform.

        This method implements the inverse FFT (i.e. from spectral space
        to real space). The values are normalised such that the spectral
        coefficient corresponds to a wave of exactly that amplitude, rather
        than the normalised value returned by the standard FFT. The sign matrix
        is used to yield a spectrum running from -N/2 to N/2 -1.

        **Parameters**

        - `array_in` : the array of spectral coefficients

        **Returns**

        - `out` : the real values at gridpoints, in the chosen framework

        **See Also**

        Note 2 in the class header (sign matrix), forward_fft

        """
        # Find side length, as spectrum may or may not have been doubled
        # for aliasing control
        side = array_in.shape[0]

        out = side**2*self.sign_mat[0:side, 0:side]*np.fft.ifft2(array_in)
        return out

    def calc_derivative(self, array_in, direction1, direction2 = False):
        """
        Computes 1st or 2nd order derivatives in x and/or y directions.

        This method implements spectral (spatial) differentiation through the use of
        the derivative matrices, stored as attributes of the SpectralToolbox class.
        The Nyquist frequency is treated specially for odd orders of multiplication.

        **Parameters**

        - `array_in` : the array of spectral coefficients to be differentiated
        - `direction1` : the first direction to be differentiated ('x' or 'y')
        - `direction2` : the optional second direction to be differentiated ('x' or 'y')

        **Returns**

        - `out` : the differentiated spectrum

        **Notes**

        Recall that the function, u, is represented in Fourier space as:

        .. math:: u(x, y, t) = \\sum_{k_{1} = -N/2}^{N/2 -1} \\sum_{k_{2} = -N/2}^{N/2 - 1} \\hat{u}(t) \\exp((2i\\pi/L)(k_{1}x + k_{2}y))

        Then the first derivative in the x-direction, as :math:`\\hat{u}(t)` is not a function of space, is simply:

        .. math:: \\frac{\\partial u(x, y, t)}{\\partial x} = \\sum_{k_{1} = -N/2}^{N/2 -1} \\sum_{k_{2} = -N/2}^{N/2 - 1} \\hat{u}(t)  (2i\\pi k_{1}/L)  \\exp((2i\\pi/L)(k_{1}x + k_{2}y))

        Thus, differentiation in spectral space is a matter of multiplying each Fourier mode by the length factor (:math:`B=2\\pi/L`) times i times the associated wavenumber.
        To reduce loops and make the implementation simpler, these modes are contained in the derivatrive matrices, which take the form:

        .. math:: \\text{deriv_mat_x1} = \\left[\\begin{array}{cccc}
                        iBk_{1} & iBk_{1} & \cdots & iBk_{1} \\\\
                       iBk_{2} & iBk_{2} & \cdots & iBk_{2}\\\\
                       \\vdots & \\vdots & \\ddots & \\vdots \\\\
                       iBk_{n} & iBk_{n} & \cdots & iBk_{n} \\end{array}\\right]

        and:

        .. math:: \\text{deriv_mat_x2} = \\left[\\begin{array}{cccc}
                        iBk_{1} & iBk_{2} & \cdots & iBk_{n} \\\\
                       iBk_{1} & iBk_{2} & \cdots & iBk_{n}\\\\
                       \\vdots & \\vdots & \\ddots & \\vdots \\\\
                       iBk_{1} & iBk_{2} & \cdots & iBk_{n} \\end{array}\\right]

        **Example**

        | ``>> A = array([[ 1.+0.j,  1.+0.j],``
        |           ``[ 1.+0.j,  1.+0.j]])``
        |
        | ``>> st.calc_derivative(A, 'x')``
        | ``>> array([[ 0.+0.j,  0.+0.j],``
        |        ``[ 0.+0.j,  0.+0.j]])``

        """
        A = array_in.copy()
        if direction1 != direction2:
            #Remove Nyquist frequency for even sample size and odd order of differentiation
            if direction1 == 'x' or direction2 == 'x':
                A[0,:] = 0.0
            if direction1 == 'y' or direction2 == 'y':
                A[:,0] = 0.0

        # Note that 'x' corresponds to the x1 direction, and 'y' to the
        # x2 direction
        # Perform first derivative in desired direction
        if direction1 == 'x':
            out = self.deriv_mat_x1*A
        elif direction1 == 'y':
            out = self.deriv_mat_x2*A

        # Perform second derivative in desired direction
        if direction2 == 'x':
            out = self.deriv_mat_x1*out
        elif direction2 == 'y':
            out = self.deriv_mat_x2*out

        return out

    def multiply_nonlinear(self, array1, array2):
        """
        Simple method for multiplying two sets of data when nonlinearities are involved
        and transformation into real space is necessary.

        This method computes the product of two quantities defined in spectral space
        using a pseudo-spectral method, i.e. the multiplication is performed in realspace.
        The spectra are padded for aliasing control and must be of the same dimension.

        **Parameters**

        - `array_1`, `array_2` : the input spectra to be multiplied

        **Returns**

        - `f3_hat` : the spectrum obtainined from the multiplication

        **See Also**

        dealias_pad, dealias_unpad

        """
        # compute grid values via FFT
        f1_vals = self.inverse_fft(self.dealias_pad(array1))
        f2_vals = self.inverse_fft(self.dealias_pad(array2))

        # multiply in space
        f3_vals = f1_vals * f2_vals

        # compute Fourier coeffs
        f3_hat = self.dealias_unpad(self.forward_fft(f3_vals))

        return f3_hat

    def solve_inverse_laplacian(self, array_in, zero_mode_val):
        """
        Solves the inverse Laplacian problem in spectral space.

        This method implements the inverse Laplacian, i.e. it finds u in:

        .. math:: \\nabla^{2}u = f

        with a constant of integration which must be specified by the user.

        **Parameters**

        - `array_in` : the known spectrum, f
        - `zero_mode_val` : the constant of integration, corresponding to the constant Fourier mode

        **Returns**

        - `array_out` : the desired quantity, u

        **See Also**

        calc_derivative

        """
        # Initialise output array, u
        array_out = np.zeros((self.N, self.N), dtype = complex)
        N = self.N

        for k1 in range(-N//2, N//2):
            for k2 in range(-N//2, N//2):
                if k1 == k2 and k1 == 0:  # Set zero-th mode, i.e. constant of integration
                    array_out[k1 + N//2, k2 + N//2] = zero_mode_val
                else:  # All other modes are found by inverting the spectral differentiation factor
                    laplace_op = -(k1**2 + k2**2)*self.factor**2
                    array_out[k1 + N//2, k2 + N//2] = array_in[k1 + N//2, k2 + N//2]/laplace_op

        return array_out

    def jacobian1(self,A):
        """
        Computes the Jacobian for a single scalar-valued function.

        This method implements the Jacobian for a single function, such as a streamfunction
        or velocity potential. The Jacobian in this case takes the form:

        .. math:: J(u) = u_{xx}u_{yy} - 2u_{xy}

        where subscripts denote partial derivatives.

        **Parameters**

        - `A` : The function of which the Jacobian is desired

        **Returns**

        - `J` : The Fourier modes of the Jacobian

        """

        # Compute second derivatives in spectral space
        A_x_x_hat = self.calc_derivative(A, 'x', 'x')
        A_y_y_hat = self.calc_derivative(A, 'y', 'y')
        A_x_y_hat = self.calc_derivative(A, 'x', 'y')
        A_y_x_hat = self.calc_derivative(A, 'y', 'x')

        # Compute realspace representations for multiplication
        A_x_x = self.inverse_fft(self.dealias_pad(A_x_x_hat))
        A_y_y = self.inverse_fft(self.dealias_pad(A_y_y_hat))
        A_x_y = self.inverse_fft(self.dealias_pad(A_x_y_hat))
        A_y_x = self.inverse_fft(self.dealias_pad(A_y_x_hat))

        # Multiply in realspace
        J_canonical = (A_x_x*A_y_y) - (A_x_y*A_y_x)

        # Return to Fourier space and return spectrum
        return self.dealias_unpad(self.forward_fft(J_canonical))

    def jacobian(self, A, B):
        """
        Computes the Jacobian for a two scalar-valued functions.

        This method implements the Jacobian for two scalar-valued functions. The Jacobian
        may be written as:

        .. math:: J(A, B) = A_{x}B_{y} - A_{y}B_{x}

        where subscripts denote partial derivatives.

        **Parameters**

        - `A`, `B` : The Fourier representations of the functions A and B (above). Must have the same dimensions.

        **Returns**

        - `J` : The Fourier modes of the Jacobian

        """

        # Compute the derivatives spectrally
        A_x_hat = self.calc_derivative(A, 'x')
        A_y_hat = self.calc_derivative(A, 'y')
        B_x_hat = self.calc_derivative(B, 'x')
        B_y_hat = self.calc_derivative(B, 'y')

        # Compute the values in realspace for multiplication
        A_x = self.inverse_fft(self.dealias_pad(A_x_hat))
        A_y = self.inverse_fft(self.dealias_pad(A_y_hat))
        B_y = self.inverse_fft(self.dealias_pad(B_y_hat))
        B_x = self.inverse_fft(self.dealias_pad(B_x_hat))

        # Compute the Jacobian
        J_canonical = (A_x*B_y) - (B_x*A_y)

        # Return to spectral space the return
        return self.dealias_unpad(self.forward_fft(J_canonical))

    def laplacian(self, array_in):
        """
        Computes the Laplacian of a given function in spectral space.

        The Laplacian may be written as:

        .. math:: f = \\nabla^{2}u = \\nabla \\cdot \\nabla u

        **Parameters**

        - `array_in` : the value of interest, in Fourier space

        **Returns**

        - `array_out` : the Laplacian of the input array

        **See Also**

        calc_derivative

        """

        # Call-through to Laplacian operator, already computed
        return self.laplace_op*array_in

#Unit Tests Below#

def testfunc(x,y):

    return np.sin(2.0*x)*np.cos(2.0*y) + np.cos(8.0*x)*np.sin(11.0*y)

def check_jacobian(st, N, L):
    """

    """
    u = np.zeros((N,N),dtype=complex)
    u_analytical = np.zeros((N,N))

    for n1 in range(0,N):
        for n2 in range(0,N):
            x = L*n1/N
            y = L*n2/N
            u[n1,n2] = testfunc(x,y)
            u_analytical[n1,n2] = (-4.0*np.sin(2.0*x)*np.cos(2.0*y) - 64.0*np.cos(8.0*x)*np.sin(11.0*y))*\
                                  (-4.0*np.sin(2.0*x)*np.cos(2.0*y) - 121.0*np.cos(8.0*x)*np.sin(11.0*y)) - \
                                  (-4.0*np.cos(2.0*x)*np.sin(2.0*y) - 88.0*np.sin(8.0*x)*np.cos(11.0*y))**2

    u_hat = st.forward_fft(u)
    u_hat_pad = st.jacobian1(u_hat)
    u_num_pad = st.inverse_fft(u_hat_pad)

    return(np.max(np.absolute(u_num_pad - u_analytical)))

def check_laplacian(st, N, L):
    """

    """
    u = np.zeros((N,N))
    u_analytical = np.zeros((N,N))

    for n1 in range(0,N):
        for n2 in range(0,N):
            x = L*n1/N
            y = L*n2/N
            u[n1,n2] = testfunc(x,y)
            u_analytical[n1,n2] = -8.0*np.sin(2.0*x)*np.cos(2.0*y) - 185.0*np.cos(8.0*x)*np.sin(11.0*y)

    u_hat = st.forward_fft(u)
    u_hat = st.laplacian(u_hat)
    u_num = st.inverse_fft(u_hat)

    return(np.max(np.absolute(u_num - u_analytical)))#/(u_analytical)

def check_derivatives(st, N, L):
    """

    """
    u = np.zeros((N,N),dtype=complex)
    u_analytical_x = np.zeros((N,N))
    u_analytical_y = np.zeros((N,N))

    for n1 in range(0,N):
        for n2 in range(0,N):
            x = L*n1/N
            y = L*n2/N
            u[n1,n2] = testfunc(x,y)
            u_analytical_x[n1,n2] = 2.0*np.cos(2.0*x)*np.cos(2.0*y) - 8.0*np.sin(8.0*x)*np.sin(11.0*y)
            u_analytical_y[n1,n2] = -2.0*np.sin(2.0*x)*np.sin(2.0*y) + 11.0*np.cos(8.0*x)*np.cos(11.0*y)

    u_hat = st.forward_fft(u)
    u_hat_x = st.calc_derivative(u_hat, 'x')
    u_hat_y = st.calc_derivative(u_hat, 'y')
    u_num_x = st.inverse_fft(u_hat_x)
    u_num_y = st.inverse_fft(u_hat_y)

    return(np.max(np.absolute(u_num_x - u_analytical_x)), np.max(np.absolute(u_num_y - u_analytical_y)))

def check_multiply(st, N, L):
    """

    """
    u1 = np.zeros((N,N),dtype=complex)
    u2 = np.zeros((N,N),dtype=complex)
    u_analytical = np.zeros((N,N))

    for n1 in range(0,N):
        for n2 in range(0,N):
            x = L*n1/N
            y = L*n2/N
            u1[n1,n2] = np.sin(2.0*x)*np.cos(4.0*y)
            u2[n1,n2] = np.sin(3.0*y)*np.cos(3.0*x)

    u_analytical = u1*u2
    u_hat1 = st.forward_fft(u1)
    u_hat2 = st.forward_fft(u2)
    u_hat = st.multiply_nonlinear(u_hat1, u_hat2)
    u_num = st.inverse_fft(u_hat)

    return(np.max(np.absolute(u_num - u_analytical)))#/(u_analytical)


if __name__ == "__main__":
    N = 256
    L = 2.0*np.pi*2.0
    st = SpectralToolbox(N, L)
    out = check_derivatives(st, N, L)
    print("Error in x-derivative is: {}".format(out[0]))
    print("Error in y-derivative is: {}".format(out[1]))
    out = check_laplacian(st, N, L)
    print("Error in Laplacian is: {}".format(out))
    out = check_jacobian(st, N, L)
    print("Error in Jacobian is: {}".format(out))
    out = check_multiply(st, N, L)
    print("Error in Multiply is: {}".format(out))

