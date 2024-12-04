# Vascular trees in Julia

Code vascular networks of the liver and the creation and simulation of 
computational models from these graphs in julia.

The corresponding python code is located in
https://github.com/matthiaskoenig/vascular_trees.

# Installation




## Install SyntheticVascularTrees.jl
This package depends on the `SyntheticVascularTrees.jl` package developed by Etienne Jessen and [Dominik Schillinger](https://www.bauing.tu-darmstadt.de/fachbereich_bau_umwelt/ueber_den_fachbereich/professoren____innen_am_fb/professorenliste_details_115392.de.jsp) at the TU Darmstadt available from 
https://gitlab.com/etiennejessen/SyntheticVascularTrees.jl.git.

- clone the repository and open terminal in `SyntheticVascularTrees.jl`
- julia --> ] --> activate . --> build --> test
- choose the Julia environment in visual studio code

# Install additional packages
```
 import Pkg; Pkg.add("OrdinaryDiffEq")
 ```


# License
* Source Code: `LGPLv3 <http://opensource.org/licenses/LGPL-3.0>`__
* Documentation: `CC BY-SA 4.0 <http://creativecommons.org/licenses/by-sa/4.0/>`__

The source is released under both the GPL and LGPL licenses version 2 or
later. You may choose which license you choose to use the software under.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License or the GNU Lesser General Public
License as published by the Free Software Foundation, either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

# Funding
Mariia Myshkina and Matthias König are supported by the Federal Ministry of Education and Research (BMBF, Germany) within the ATLAS project (grant number 031L0304B).

Matthias König is supported by the Federal Ministry of Education and Research (BMBF, Germany)
within the research network Systems Medicine of the Liver (**LiSyM**, grant number 031L0054) 
and by the German Research Foundation (DFG) within the Research Unit Programme FOR 5151 
"`QuaLiPerF <https://qualiperf.de>`__ (Quantifying Liver Perfusion-Function Relationship in Complex Resection - A Systems Medicine Approach)" by grant number 436883643 and by grant number 465194077 (Priority Programme SPP 2311, Subproject SimLivA).


© 2024 Mariia Myshkina & Matthias König