# -*- python -*-
# $Header$
# Authors: Tsunefumi Mizuno <suhonen@slac.stanford.edu>
# Version: CRflux-01-17-05
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('addLinkDeps', package='CRflux', toBuild='component')
CRflux=libEnv.SharedLibrary('CRflux', listFiles(['src/*.cxx','src/psb97/*.cxx',
                                                 'src/Dll/*.cxx']))

progEnv.Tool('CRfluxLib')

test_CRflux = progEnv.GaudiProgram('test_CRflux',listFiles(['src/test/*.cxx']),
                                   test = 1, package='CRflux')

#if baseEnv['PLATFORM'] != 'win32':
progEnv.Tool('registerTargets', package = 'CRflux',
             libraryCxts = [[CRflux, libEnv]],
             testAppCxts = [[test_CRflux, progEnv]],
             includes = listFiles(['CRflux/*.h']),
             jo=['src/test/jobOptions.txt'])
	





