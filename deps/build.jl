using BinDeps

import BinDeps.library_dependency

@BinDeps.setup

deps = [
		metis = library_dependency("metis",aliases = ["libmetis"])
		dmumps = library_dependency("libdmumps")
		zmumps = library_dependency("libzmumps")
		commonmumps = library_dependency("libmumps_common")
	   ]

# download sources
provides(Sources,
        {
		  URI("http://mumps.enseeiht.fr/MUMPS_4.10.0.tar.gz") => [dmumps,zmumps,commonmumps],
		  URI("http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.3.tar.gz") => metis
		})

# build from source
provides(BuildProcess, 
			{ 
				Autotools(libtarget = "") => metis,
				Autotools(libtarget = "d") => [dmumps,commonmumps],
				Autotools(libtarget = "z") =>zmumps
			})

@BinDeps.install