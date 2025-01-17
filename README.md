[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Incorporating promotional effects in sales planning of retail industry using geometric programming

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The data in this repository was used in the research reported on in the paper
[Incorporating promotional effects in sales planning of retail industry using geometric programming](https://doi.org/10.1287/ijoc.2023.0275) by M. Khandan, and P. Hoseinpour.

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2023.0275

https://doi.org/10.1287/ijoc.2023.0275.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{Khandan2024ijoc,
  author =        {Khandan, Melika and Hoseinpour, Pooya},
  publisher =     {INFORMS Journal on Computing},
  title =         {{Incorporating promotional effects in sales planning of retail industry using geometric programming}},
  year =          {2024},
  doi =           {10.1287/ijoc.2023.0275.cd},
  url =           {https://github.com/INFORMSJoC/2023.0275},
  note =          {Available for download at https://github.com/INFORMSJoC/2023.0275},
}
```

## Description

This directory contains all the instances and results in both **Numerical Experiments** and **Case Study** Sections in the paper. The directory contains the folders `Instances for Numerical Experiments (Section 5)` and `Instances for Case study (Section 6)`:
- `Instances for Numerical Experiments (Section 5)`: includes the instances and results reported in Section 5. This folder is organized as follows:
  - `Instances for Numerical Experiments (Section 5)/for Comparison between convex and non-convex model (Section 5-1)/`: instances and results reported in Table 1.
  - `Instances for Numerical Experiments (Section 5)/for Lagrangian Decomposition Algorithm (Section 5-2)/`: instances and results reported in Section 5-2.
    - `Instances for Numerical Experiments (Section 5)/for Lagrangian Decomposition Algorithm (Section 5-2)/Table 5/`: instances and results reported in Table 5.
    - `Instances for Numerical Experiments (Section 5)/for Lagrangian Decomposition Algorithm (Section 5-2)/Table 6 and Table 7/`: instances and results reported in Table 6 and Table 7.
  - `Instances for Numerical Experiments (Section 5)/for Sensitivity Analysis (Section 5-3)/`: instances and results reported in Section 5-3.
    - `Instances for Numerical Experiments (Section 5)/for Sensitivity Analysis (Section 5-3)/Table 9 and Figure 5/`: instances and results reported in Table 9 and Figure 5.
- `Instances for Case study (Section 6)`: includes the instances and results reported in Section 6. This folder is organized as follows:
  - `Instances for Case study (Section 6)/Figure 6 and Table 11/`: instances and results reported in Figure 6 and Table 11.
  - `Instances for Case study (Section 6)/Figure 7/`: instances and results reported in Figure 7.
  - `Instances for Case study (Section 6)/Figure 8/`: instances and results reported in Figure 8.

In addition, the **source code** is availabe, named as "main code.gms", implemented in GAMS. The file includes the following four parts:
1. Codes for the convex promotion optimization model using the geometric programming approach.
2. Codes for the non-convex promotion optimization model. This part is added to ensure the accuracy of the Taylor approximation based on the comments we received in the second revision.
3. Codes for the Lagrangian decomposition algorithm.
4. Codes for the Lagrangian decomposition algorithm, decomposed over time (t) to reduce the running time of sub-model 1.
