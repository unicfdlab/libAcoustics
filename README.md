
Critical parameters for FfowcsWilliams-Hawkings analogy.

The proper choice is the CELL interpolation scheme. Other ones will lead to errors.

interpolationScheme. choice of
      cell          : use cell-centre value only; constant over cells (default)
      cellPoint     : use cell-centre and vertex values
      cellPointFace : use cell-centre, vertex and face values.
      pointMVC      : use point values only (Mean Value Coordinates)
      cellPatchConstrained : like 'cell' but uses cell-centre except on
                             boundary faces where it uses the boundary value.
                             For use with e.g. patchCloudSet.
1] vertex values determined from neighbouring cell-centre values
2] face values determined using the current face interpolation scheme for the field (linear, gamma, etc.)