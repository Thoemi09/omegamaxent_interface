.. _UserGuide:

User Guide
==========

This is the user guide for the OmegaMaxEnt_TRIQS module. It should be sufficient ot obtain good results in many cases. For a more advanced use, which will sometimes be necessary, refer to the :math:`\Omega MaxEnt` `user guide`_.

Main function in OmegaMaxEnt_TRIQS
----------------------------------

Whether you have a scalar-, a matrix-valued or a block Green's function, the real frequency Green's function is obtained with the function **compute_GfReFreq()**, which takes an object of type Gf,GfImFreq, GfImTime, or a BlockGf containing objects of one of those types, and returns an object of type GfReFreq or a BlockGf containing GfReFreq objects. The function signature is::

    compute_GfReFreq(G, ERR=None, grid_params=[], name="$G^R(\omega)$", interactive_mode=True,
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

grid_params:
    Optional list of the form :math:`[\omega_{min}, \Delta\omega, \omega_{max}]`.

    Defines the minimum frequency, the frequency step, and the maximum frequency, respectively, of the real frequency grid of the output Green function.

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
    Optional list of the form [:math:`\Delta\omega`] or [:math:`\Delta\omega`, spectrum_width]` or [:math:`\Delta\omega`, spectrum_width, spectrum_center].

    Grid parameters used in the computation.

    :math:`\Delta\omega` is the frequency step used in the main spectral region, namely, the part of the grid where most of the spectral weight is located.

    spectrum_width is the width of the main spectral region (should be between 2 and 4 standard deviations typically).

    spectrum_center is the center of the main spectral region.

non_uniform_grid:
    Optional boolean.

    Tells :math:`\Omega MaxEnt` to use a non-uniform grid in the main spectral region for the computation.

    This accelerates the calculation if the spectrum has a peak at zero frequency that is very narrow compared to the total width of the spectrum.

inv_sym:
    Optional boolean.

    If G is a matrix or a BlockGf, set ``inv_sym=True`` if :math:`G_{ji}=G_{ij}`. This simplifies the calculation of the off diagonal elements.

mu, nu:
    Optional parameters involved in the calculation of off-diagonal elements of matrix-valued Green functions. See appendix C of the :math:`\Omega MaxEnt` `user guide`_ for more details.


:math:`\Omega MaxEnt` parameter files
--------------------------------------

:math:`\Omega MaxEnt` uses the file **OmegaMaxEnt_input_params.dat** to interact with the user. You can optionally create that file in advance with the function **create_params_file()**. This allows you to set some parameters that are not set by **compute_GfReFreq()**. Some parameters, which appear in the top section of the file, are exclusively set by **compute_GfReFreq()**, but all the other parameters can be modified.

If **OmegaMaxEnt_input_params.dat** does not exist when **compute_GfReFreq()** is called, it will create it.

:math:`\Omega MaxEnt` also uses a file called **OmegaMaxEnt_other_params.dat**, also created by **create_params_file()**, which defines a certain number of internal parameters on which the computation depends. Do not modify this file unless you are an advanced user.

All the parameters in **OmegaMaxEnt_input_params.dat** and **OmegaMaxEnt_other_params.dat** are described in the `user guide`_.

.. _`interactive mode`:


Interactive mode
----------------

In interactive mode, :math:`\Omega MaxEnt` displays figures during the execution. If parameter "display preprocessing figures" in **OmegaMaxEnt_input_params.dat** is enabled, figures are displayed during the preprocessing stage. Otherwise, figures are displayed only at the end of the calculation, showing the resulting Green function, along with different quantities used as diagnostic tools. Using those tools is very useful to assess, first, if the result is valid and, second, if it is the best result possible given the data. Therefore, when processing a set of data for the first time, it is strongly advised to use the interactive mode. Details about how to interpret the diagnostic quantities are given in the `user guide`_.

In interactive mode, :math:`\Omega MaxEnt` pauses at the end of the preprocessing stage if parameter "preprocess only" in **OmegaMaxEnt_input_params.dat** is enabled, and at the end of the calculation. During a pause, you can modify **OmegaMaxEnt_input_params.dat** and resume the execution of :math:`\Omega MaxEnt` once the file is saved. Otherwise, if the calculation is over and you are satisfied with the result displayed, you can exit the execution by closing all the figures and entering any character other than ``'y'`` in the terminal. This will resume the execution of the python function **compute_GfReFreq()**.

If ``interactive_mode=False``, :math:`\Omega MaxEnt` will not display any figure and will exit at the end of the calculation, resuming the execution of **compute_GfReFreq()**.

.. note::

    For the continuation of **matrix-valued** or **block** Green's functions, :math:`\Omega MaxEnt` is called  the same number of times as there are elements in each matrix (or in the upper part if ``inv_sym=True``). If you are in interactive mode, figures showing the result will appear each time and, once you have closed them, you have to tell the program **not** to continue execution to let the analytic continuation of the matrix or block function continue.


Imaginary time data
-------------------

If your data is a scalar GfImTime and you do not have an estimate of the error, or the error is constant, do not set parameter ``ERR``. Otherwise, because :math:`\Omega MaxEnt` works internally in Matsubara frequency, it will Fourier transform the covariance matrix, which is not useful in that case because the result will also be a constant diagonal covariance in frequency and the result does not depend on the absolute value of the error. Avoiding the Fourier transform of the covariance matrix will therefore save computation time without changing the result.

On the other hand, if the error depends on :math:`\tau` and you use ``ERR`` to provide it, note that the Fourier transform of the Green function is saved by default as a GfImFreq object called 'G' in file "G_im_freq.h5" and the Fourier transform of the covariance matrix is saved in files "covar_ReRe.dat", "covar_ImIm.dat" and "covar_ReIm.dat" in directory "Fourier_transformed_data". This can be useful if you want to perform the continuation again on the same data. Then you can pass the saved GfImFreq object to **compute_GfReFreq()** instead of the original GfImTime object and use the parameters "re-re covariance file", "im-im covariance file" and "re-im covariance file" in section INPUT FILES PARAMETERS of the file **OmegaMaxEnt_input_params.dat** to provide the covariance to :math:`\Omega MaxEnt`.


Output Figures
--------------

If ``save_figures_data=True``, reagrdless of the value of ``interactive_mode``, you can display the same figures that are displayed in interactive mode with the function **display_figures()** after the execution of **compute_GfReFreq()**. Note however that only the figures for the last continuation done by :math:`\Omega MaxEnt` in a given directory are accessible.


Choice of frequency grid
------------------------




For more details on how to use :math:`\Omega MaxEnt`, see the `user guide`_.


Example: Suppose you have saved a TRIQS Matsubara Green's function as 'G' in file "G.h5", here is a script to obtain the corresponding real frequency Green's function::


    from pytriqs.archive import HDFArchive as HA
    import OmegaMaxEnt_TRIQS as OT

    #load the Green's function
    A=HA("G.h5",'r'):
    G=A['G']

    GR=OT.compute_GfReFreq(G)



.. _`user guide`: https://www.physique.usherbrooke.ca/MaxEnt/index.php/User_Guide