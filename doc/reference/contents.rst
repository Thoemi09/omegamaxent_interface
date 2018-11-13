.. _UserGuide:

User Guide
==========

This is the user guide for the OmegaMaxEnt_TRIQS module. It should be sufficient ot obtain good results in many cases. For a more advanced use, which may sometimes be required, refer to the :math:`\Omega MaxEnt` `user guide`_.


In the following, we refer to the function to be analytically continued as the "Green's function", but it can be also a two-particle correlation function, a self-energy, or any function having a representation of the form

.. math::

    \mathbf{G}(i\omega_n)=\int_{-\infty}^{\infty} d\omega \frac{\mathbf{A}(\omega)}{i\omega_n-\omega},

where :math:`\omega_n` is a fermionic or boson Matsubara frequency, or

.. math::

    \mathbf{G}(\tau)=\int_{-\infty}^{\infty} d\omega \frac{e^{-\omega\tau}\mathbf{A}(\omega)}{1\pm e^{-\beta\omega}},

where the + and - signs apply to the fermionic and bosonic cases, respectively, and where :math:`\mathbf{A}(\omega)` is a positive semi-definite matrix in the fermionic case or :math:`\mathbf{A}(\omega)/\omega` is positive semi-definite in the bosonic case.

Obtain the retarded Green's function
------------------------------------

Whether you have a scalar-, a matrix-valued or a block Green's function, the real frequency Green's function is obtained with the function **compute_GfReFreq()**, which takes an object of type Gf, GfImFreq, GfImTime, or a BlockGf containing objects of one of those types, and returns an object of type GfReFreq or a BlockGf of GfReFreq objects. The function signature is::

    compute_GfReFreq(G, ERR=None, grid_params=[], name="$G^R$", interactive_mode=True,
                     save_figures_data=True, save_G=True, comp_grid_params=[], non_uniform_grid=False,
                     inv_sym=False, mu=1, nu=1)

where the parameters are:

G:
    Gf, GfImFreq, GfImTime or BlockGf object.

    The input Matsubara Green function.

ERR:
    Optional numpy array.

    Standard deviation if G is scalar or has a single element.

    ERR must have the same shape as G.data.

    For a non-diagonal covariance, see the :math:`\Omega MaxEnt` `user guide`_.

.. _grid_params:

grid_params:
    Optional list of the form :math:`[\omega_{min}, \Delta\omega, \omega_{max}]`.

    Defines the minimum frequency, the frequency step, and the maximum frequency, respectively, of the frequency grid on which the output Green function will be defined.

    If empty, the output grid is set by :math:`\Omega MaxEnt`.

    Note: this is only the output frequency grid. For the grid used in the calculation, see parameter comp_grid_params_.

name:
    Optional string.

    Name parameter of the returned GfReFreq object

interactive_mode:
    Optional boolean.

    Turns off the `interactive mode`_ of :math:`\Omega MaxEnt` if set to False.

save_figures_data:
    Optional boolean.

    Tells :math:`\Omega MaxEnt` not to save the figures files if set to False.

    Set to False if G is a matrix or a BlockGf.

save_G:
    Optional boolean.

    By default, the result is save in hdf5 format in file "G_Re_Freq.h5".

.. _comp_grid_params:

comp_grid_params:
    Optional list of the form [:math:`\Delta\omega`] or [:math:`\Delta\omega`, SW] or [:math:`\Delta\omega`, SW, SC].

    Grid parameters used in the computation.

    :math:`\Delta\omega` is the frequency step used in the main spectral region, namely, the part of the grid where most of the spectral weight is located.

    SW is the width of the main spectral region (should be between 2 and 4 standard deviations typically).

    SC is the center of the main spectral region.

non_uniform_grid:
    Optional boolean.

    Tells :math:`\Omega MaxEnt` to use a non-uniform grid in the main spectral region for the computation.

    This accelerates the calculation if the spectrum has a peak at zero frequency that is very narrow compared to the total width of the spectrum.

inv_sym:
    Optional boolean.

    If G is a matrix or a BlockGf, set ``inv_sym=True`` if :math:`G_{ji}=G_{ij}`. This simplifies the calculation of the off diagonal elements.

.. _`mu and nu`:

mu, nu:
    Optional parameters involved in the calculation of off-diagonal elements of matrix-valued Green functions. See appendix C of the :math:`\Omega MaxEnt` `user guide`_ for more details.


:math:`\Omega MaxEnt` parameter files
--------------------------------------

:math:`\Omega MaxEnt` uses the file **OmegaMaxEnt_input_params.dat** to interact with the user. You can optionally create that file in advance with the function **create_params_file()**. This allows you to set some parameters that are not set by **compute_GfReFreq()**. Some parameters, which appear in the top section of the file, are exclusively set by **compute_GfReFreq()**, but all the other parameters can be modified.

If **OmegaMaxEnt_input_params.dat** does not exist when **compute_GfReFreq()** is called, it will create it.

:math:`\Omega MaxEnt` also uses a file called **OmegaMaxEnt_other_params.dat**, also created by **create_params_file()**, which defines a certain number of internal parameters on which the computation depends. Do not modify this file unless you are an advanced user.

All the parameters in **OmegaMaxEnt_input_params.dat** and **OmegaMaxEnt_other_params.dat** are described in the :math:`\Omega MaxEnt` `user guide`_.

.. _`interactive mode`:


Interactive mode
----------------

In interactive mode, :math:`\Omega MaxEnt` displays figures during the execution. If parameter "display preprocessing figures" in **OmegaMaxEnt_input_params.dat** is enabled, figures are displayed during the preprocessing stage. Otherwise, figures are displayed at the end of the calculation, showing the resulting Green function, along with different quantities used as diagnostic tools. Using those tools is very useful to assess, first, if the result is valid and, second, if it is the best result possible given the data. Therefore, when processing a set of data for the first time, it is strongly advised to use the interactive mode. Details about how to interpret the diagnostic quantities are given in the :math:`\Omega MaxEnt` `user guide`_.

You can also force the calculation to pause using parameters "minimum value of alpha" or "number of values of alpha computed in one execution" in **OmegaMaxEnt_input_params.dat** and look at the results at different stages (i.e. values of alpha_). You can always modify **OmegaMaxEnt_input_params.dat** during those pauses and the changes will be applied when execution is resumed at the point of interuption.

In interactive mode, :math:`\Omega MaxEnt` pauses at the end of the preprocessing stage if parameter "preprocess only" in **OmegaMaxEnt_input_params.dat** is enabled, and at the end of the calculation. During a pause, you can modify **OmegaMaxEnt_input_params.dat** and resume the execution of :math:`\Omega MaxEnt` once the file is saved. Otherwise, if the calculation is over and you are satisfied with the result displayed, you can exit the execution by closing all the figures and entering any character other than ``'y'`` in the terminal. This will resume the execution of the python function **compute_GfReFreq()**.

If ``interactive_mode=False``, :math:`\Omega MaxEnt` will not display any figure and will exit at the end of the calculation, resuming the execution of **compute_GfReFreq()**.

If you find that there are too many figures, instead of completely disable the interactive mode to eliminate all of them, you can reduce the number of figures displayed by turning off only some of them in section DISPLAY OPTIONS of **OmegaMaxEnt_input_params.dat**.

.. note::

    For the continuation of **matrix-valued** Green's functions, :math:`\Omega MaxEnt` is called  the same number of times as there are elements in the matrix (or in the upper part if ``inv_sym=True``). If you are in interactive mode, figures showing the result will appear each time and, once you have closed them, you have to tell the program **not** to continue execution to let the analytic continuation of the matrix or block function continue.


Imaginary time data
-------------------

If your data is a scalar GfImTime and you do not have an estimate of the error, or the error is constant, do not set parameter ``ERR``. Otherwise, because :math:`\Omega MaxEnt` works internally in Matsubara frequency, it will Fourier transform the covariance matrix, which is not useful in that case because the result will also be a constant diagonal covariance in frequency, and the result does not depend on the absolute value of the error. Avoiding the Fourier transform of the covariance matrix will therefore save computation time without changing the result.

On the other hand, if the error depends on :math:`\tau` and you use ``ERR`` to provide it, note that the Fourier transform of the Green function is saved by default as a GfImFreq object called 'G' in file "G_im_freq.h5" and the Fourier transform of the covariance matrix is saved in files "covar_ReRe.dat", "covar_ImIm.dat" and "covar_ReIm.dat" in directory "Fourier_transformed_data". This can be useful if you want to perform the continuation again on the same data. Then you can pass the saved GfImFreq object to **compute_GfReFreq()** instead of the original GfImTime object and use the parameters "re-re covariance file", "im-im covariance file" and "re-im covariance file" in section INPUT FILES PARAMETERS of the file **OmegaMaxEnt_input_params.dat** to provide the covariance to :math:`\Omega MaxEnt`.


Display figures after execution
-------------------------------

If ``save_figures_data=True``, regardless of the value of ``interactive_mode``, you can display the same figures that are displayed in interactive mode with the function **display_figures()** after the execution of **compute_GfReFreq()**. Note however that only the figures for the last continuation done by :math:`\Omega MaxEnt` in a given directory are accessible. More details on the output figures are given in the :math:`\Omega MaxEnt` `user guide`_.


.. _real_freq_grid:

Frequency grids
----------------

There are two parameters related to the real frequency grid: grid_params_ and comp_grid_params_.

Parameter grid_params_, allows you to control the frequency grid of the output Green's function, which is a uniform density grid defined between :math:`\omega_{min}` and :math:`\omega_{max}` with step :math:`\Delta\omega`. This is an optional parameter. If not provided, :math:`\Omega MaxEnt` generates an output frequency grid that is usually well adapted to the spectrum.

For computational efficiency reasons, the real frequency grid used in the calculation is different from the output grid. In many cases the default computational grid generated by :math:`\Omega MaxEnt` is well adapted to the spectrum and there is no need for the user to set parameter comp_grid_params_. In case the result is not satisfactory, you can use that parameter to control the computational grid in the region where most of the spectral weight is located. Outside that region, a particular non-uniform grid is always used by :math:`\Omega MaxEnt`. More advanced parameters are also available to control the computational grid in section FREQUENCY GRID PARAMETERS of file **OmegaMaxEnt_input_params.dat**. See the :math:`\Omega MaxEnt` `user guide`_ for more details on those parameters.

.. _alpha:

Choice of entropy weight :math:`\alpha`
---------------------------------------

In the `maximum entropy`_ method, there is an adjustable parameter :math:`\alpha`, i.e. the weight of the entropy term, which can be chosen in `different ways`_. :math:`\Omega MaxEnt` chooses the value where the curvature of :math:`log(\chi^2)` as a function of :math:`\gamma log(\alpha)` is maximal [#OME]_. Here :math:`\gamma<1` is a parameter that reduces the probability of a wrong value of :math:`\alpha` to be chosen (:math:`\gamma` is defined in **OmegaMaxEnt_other_params.dat**, default value: :math:`\gamma=0.2`). Despite the use of :math:`\gamma` and some smoothing of the curve :math:`log(\chi^2)` vs :math:`\gamma log(\alpha)` in the computation of the curvature, there is still a chance that a wrong value of :math:`\alpha` will be selected because of some irregularities in :math:`log(\chi^2)` vs :math:`\gamma log(\alpha)` that produce parasitic peaks in the curvature. This is one of the reasons why the diagnostic tools are useful.

Matrix-valued functions
------------------------------

If the Green's function is matrix-valued, the calculation is done using the auxiliary Green's function approach described in [#AuxME]_ or appendix C of the :math:`\Omega MaxEnt` `user guide`_. In that calculation, the off-diagonal elements of the retarded function are obtained indirectly from the spectral functions of the diagonal elements and auxiliary functions. Those auxiliary Matsubara functions are constructed to have positive semi-definite spectral functions, so that they can be computed with the standard maximum entropy approach. As for the diagonal elements, they always have positive semi-definite spectral functions. Therefore, in that calculation, the retarded functions corresponding to the diagonal and the auxiliary Matsubara functions are first computed with :math:`\Omega MaxEnt` and then combined at the python level to obtain the retarded off-diagonal elements.

For Green or correlation functions of the form :math:`\langle T_{\tau} o_i(\tau) o_j^\dagger\rangle`, where :math:`o_i` and :math:`o_j` are the same type of excitations, e.g. electronic excitations, the parameters `mu and nu`_ should be left equal to 1. On the other hand, for a correlation function of the form :math:`\langle T_{\tau} p(\tau) q^\dagger\rangle`, where :math:`p` and :math:`q` are different types of excitations, different values of `mu and nu`_ should be tried to find a stable result [#AuxME]_.


Example
-------

Suppose you have saved a Matsubara Green's function as a TRIQS object 'G' in a hdf5 file "G.h5". Here is the simplest way to obtain the corresponding real frequency Green's function::

    from pytriqs.archive import HDFArchive as HA
    import OmegaMaxEnt_TRIQS as OT

    #load the Matsubara Green's function
    with HA("G.h5",'r') as A:
        G=A['G']

    #obtain the retarded Green's function
    GR=OT.compute_GfReFreq(G)

.. _`user guide`: https://www.physique.usherbrooke.ca/MaxEnt/index.php/User_Guide
.. _`maximum entropy`: https://triqs.github.io/maxent/jenkins/basicnotions/mathematics.html
.. _`different ways`: https://triqs.github.io/maxent/jenkins/basicnotions/maxentflavors.html#maxent-flavors

.. [#OME] `D.Bergeron and A.-M.S. Tremblay. Phys. Rev. E, 94:023303, 2016 <https://journals.aps.org/pre/abstract/10.1103/PhysRevE.94.023303>`_

.. [#AuxME] `A. Reymbaut, A.-M. Gagnon, D. Bergeron, A.-M. S. Tremblay. Phys. Rev. B 95:121104, 2017 <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.95.121104>`_
