# CppBOLOS

CppBOLOS is an open-source C++ version of the original [BOLOS](https://github.com/aluque/bolos/tree/master#hp2005), 
an electron Boltzmann equation solver.

## Overview

The original [BOLOS](https://github.com/aluque/bolos/tree/master#hp2005) library is a pure Python solution for the Boltzmann equation for electrons in a non-thermal plasma. 
It builds on prior work, primarily by G. J. M. Hagelaar and L. C. Pitchford, who developed [Bolsig+](https://www.bolsig.laplace.univ-tlse.fr/index.html). 
The CppBOLOS project maintains compatibility with the Bolsig+ cross-section input format and offers enhanced performance benefits due to its C++ implementation.
This was achieved with the assistance of OpenAI's [GPT-4](https://openai.com/gpt-4) model, which played a crucial role in the translation process.
CppBOLOS is suitable for integration into various other C++ scientific codes, such as Cantera and OpenFOAM. This facilitates comprehensive studies of multi-physics problems involving low-temperature plasma.

The motivation behind this translation was:

1. Leverage the speed and efficiency of C++ while maintaining the core functionalities and logic of the original BOLOS (two-term expansion).
2. Explore the potential of GPT-4 for the translation/development of complex scientific code.

We have incorporated CppBOLOS into [Cantera](https://cantera.org) at source level to build [ChemPlasKin](https://doi.org/10.1016/j.jaecs.2024.100280), 
a general purpose freeware for zero-dimensional (0D) simulations of neutral gas chemical kinetics coupled with non-equilibrium plasma.

## Key Features

- **High Performance**: Leverage the speed and optimization of C++; fivefold speedup compared to BOLOS; around fourfold speedup compare to Bolsig+ (Windows version).
- **Compatibility**: Maintains compatibility with Bolsig+ cross-section input format.
- **Open Source**: Freely available for modification and use.

## Acknowledgments

- **Original BOLOS**: This project owes its foundation to the original BOLOS developed by Alejandro Luque at the Instituto de Astrofísica de Andalucía (IAA), CSIC. The original BOLOS is licensed under the LGPLv2 License and can be accessed on [GitHub](https://github.com/aluque/bolos/tree/master) and its documentation on [ReadTheDocs](http://bolos.readthedocs.org).

- **GPT-4 Assistance**: The translation and development of CppBOLOS were significantly expedited with the assistance of OpenAI's [GPT-4](https://openai.com/gpt-4) model, showcasing the capability of AI in scientific code translation.

## Getting Started

1. **Prerequisites**: 
- A C++ compiler supporting C++17 (e.g., GCC, Clang). 
- CMake (version >= 3.17).
2. **Clone the Repository**:
   
   ```sh
   git clone --recurse-submodules https://github.com/ShaoX96/CppBOLOS.git
   ```

3. **Build the Project**:
   
   ```sh
   mkdir build
   cd build
   cmake ..
   make
   ```
   
   The cross section data is stored in the `data/` directory.

## Contributing

We welcome contributions to CppBOLOS! If you'd like to contribute, please fork the repository and make changes as you'd like. Pull requests are warmly welcomed.

## License

CppBOLOS is open-sourced under the [LGPLv2 License](https://www.gnu.org/licenses/old-licenses/lgpl-2.0.html).

## Links

- [Original BOLOS GitHub](https://github.com/aluque/bolos/tree/master)
- **Bolsig+ Reference Paper**: [G. J. M. Hagelaar and L. C. Pitchford, "Solving the Boltzmann equation to obtain electron transport coefficients and rate coefficients for fluid models", Plasma Sources Science and Technology, 2005](https://iopscience.iop.org/article/10.1088/0963-0252/14/4/011)
