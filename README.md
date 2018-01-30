# CRNSimulator Mathematica Package

A Mathematica package for working with networks of coupled chemical
reactions. It forms a foundation for syntactic manipulation of
chemical reaction networks as Mathematica expressions. It allows
mixing and matching mass-action kinetics with other kinds of dynamics,
and provides a simple way to simulate experiments in which a sequence
of chemical additions is performed. It is particularly well-suited for
engineered chemical systems, in which chemical reaction networks can
be used as a kind of "programming language".

See Examples.nb for the kind of things it can do.

- Package: [CRNSimulator.m](src/core/CRNSimulator.m)
- Extensions: [CRNSimulatorExtensions.m](src/core/CRNSimulatorExtensions.m) (needed for some more advanced features)
- Examples: [Examples.nb](src/core/Examples.nb)

## CRNSimulatorSSA Mathematica Package [Beta]

Extends CRNSimulator to perform Gillespie SSA simulation. The main
simulation loop is compiled to C code for speed. Input reactions and
initial counts using the same syntax as CRNSimulator.

See ExamplesSSA.nb for simple examples.

- Package: [CRNSimulatorSSA.m](src/ssa/CRNSimulatorSSA.m) (requires CRNSimulator package)
- Examples: [ExamplesSSA.nb](src/ssa/ExamplesSSA.nb)