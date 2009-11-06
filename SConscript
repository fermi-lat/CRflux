# -*- python -*-
# $Header$
# Authors: Tsunefumi Mizuno <suhonen@slac.stanford.edu>
# Version: CRflux-01-17-02
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('CRfluxLib', depsOnly = 1)
CRflux=libEnv.SharedLibrary('CRflux', listFiles(['src/*.cxx','src/psb97/*.cxx',
                                                 'src/Dll/*.cxx']))

progEnv.Tool('CRfluxLib')
test_CRflux = progEnv.GaudiProgram('test_CRflux',listFiles(['src/test/*.cxx']),
                                   test = 1)

#if baseEnv['PLATFORM'] != 'win32':
progEnv.Tool('registerTargets', package = 'CRflux',
             libraryCxts = [[CRflux, libEnv]],
             testAppCxts = [[test_CRflux, progEnv]],
             includes = listFiles(['CRflux/*.h']))
	





