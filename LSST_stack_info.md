##Details on working with the LSST stack:
- Download here: https://community.lsst.org/t/developing-software-with-the-dm-stack/126
- To run, use these commands in the lsstsw directory:

```
$ source bin/setup.sh          (find this command using: $ ./bin/deploy)
$ setup -r ./stack/DarwinX86/afw/(*version number inserted here*)/      
```
(go to ./stack/DarwinX86/afw to see what version numbers are available)

To edit the lsst code, edit the copy in (*path to lsstsw*)/build
	To run, if for example code in (*path to lsstsw*)/build/afw has been edited:

	$ cd (*path to lsstsw*)/build/afw
	$ setup -r . -t (*appropriate version number here*)
	$ scons opt=3

	The version number should match the version of an installed copy of afw located in (*path to lsstsw*)/stack that corresponds to the state of afw in (*path to lsstsw*)/build/afw before the edits were made.  To find the versions of afw that have been installed (with example output):

	$ eups list afw
10.1-40-gcb6923d 	b1638
  	11.0-7-g974e250 	b1640   	//if this is the copy I want, I would run:
	$ cd (*path to lsstsw*)/build/afw
	$ setup -r . -t b1640   //this is where the b1640 goes
	$ scons opt=3

To call a function you created (for example terraFuncNameInC()) whose declaration is in a .h file you created (blurExampleTerraStandalone.h) with an associated .o file (blurExampleTerraStandalone.o) from an afw file (lsstsw/build/afw/src/math/detail/BasicConvolve.cc):

1.  place blurExampleTerraStandalone.h in lsstsw/build/afw/include/lsst/afw/math/Terra/blurExampleTerraStandalone.h
2.  rename blurExampleTerraStandalone.o as blurExampleTerraStandalone.os and place it in /lsstsw/build/afw/src/math/blurExampleTerraStandalone.os
	3.  add:
 #include "lsst/afw/math/Terra/blurExampleTerraStandalone.h" 
to the file lsstsw/build/afw/src/math/detail/BasicConvolve.cc (and desired call)
4.  add:
objs.append("#src/math/blurExampleTerraStandalone.os")
	to the file:
	lsstsw/build/afw/lib/SConscript
	5.  change:
scripts.BasicSConstruct("afw")
to:
env = scripts.BasicSConstruct("afw")
	env['STATIC_AND_SHARED_OBJECTS_ARE_THE_SAME']=1
	in the file:
	lsstsw/build/afw/SConstruct
	6.  As above:
	$ cd (*path to lsstsw*)/build/afw
	$ setup -r . -t (*appropriate version number here*)
	$ scons opt=3

	-This may not be the best method, but seems to work.



Running image difference pipeline on Lightroast:
Setup:
git clone https://github.com/LSST/pipe_tasks
git clone https://github.com/lsst/obs_decam

Run:
jkuck@lightroast:~/lsstsw$ ./bin/deploy

Done. Run the following:

    . /home/jkuck/lsstsw/bin/setup.sh

to begin using it.

jkuck@lightroast:~/lsstsw$ . /home/jkuck/lsstsw/bin/setup.sh
notice: lsstsw tools have been set up.
jkuck@lightroast:~/lsstsw$ eups list afw
   10.1-33-0  	b1589
   10.1-33-g17a53c7 	b1588
   10.1-33-g17a53c79 
   11.0-8-g38426eb 	b1590
jkuck@lightroast:~/lsstsw$ cd pipe_tasks/
jkuck@lightroast:~/lsstsw/pipe_tasks$ eups list pipe_tasks
   10.1-26-g073f8e3+3 	b1588
   10.1-27-gd853fd0 	b1589
   11.0-14-ga314014 	b1590
jkuck@lightroast:~/lsstsw/pipe_tasks$ setup pipe_tasks -t b1590
jkuck@lightroast:~/lsstsw/pipe_tasks$ setup -r . -k		// -r setup specified by the following path
// -k keep everything else as is except what's //specified by this path and dependencies

jkuck@lightroast:~/lsstsw/pipe_tasks$ cd ../obs_decam/
jkuck@lightroast:~/lsstsw/obs_decam$ setup -r . -k
jkuck@lightroast:~/lsstsw/obs_decam$ /home/jkuck/lsstsw/pipe_tasks/bin/imageDifference.py /home/jkuck/lsstsw/diffim_test_data  --id visit=289820 ccdnum=10 --templateId visit=288976 ccdnum=10 --configfile /home/jkuck/lsstsw/diffim_test_data/diffimconfig.py --output /home/jkuck/lsstsw/diffim_test_data/differencingOutput --clobber-config

Extra details
These python files determine pipeline setup:
/pipe_tasks/python/lsst/pipe/tasks/imageDifference.py
/lsstsw/stack/DarwinX86/ip_diffim/11.0-4-gf67ab7e/python/lsst/ip/diffim/imagePsfMatch.py

