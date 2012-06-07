# -*- python -*-
# $Header$
# Authors: Tsunefumi Mizuno <suhonen@slac.stanford.edu>
# Version: CRflux-01-18-02

Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('addLinkDeps', package='CRflux', toBuild='component')
CRflux=libEnv.ComponentLibrary('CRflux',
                               listFiles(['src/*.cxx','src/psb97/*.cxx']))

progEnv.Tool('CRfluxLib')

test_CRflux = progEnv.GaudiProgram('test_CRflux',listFiles(['src/test/*.cxx']),
                                   test = 1, package='CRflux')

#if baseEnv['PLATFORM'] != 'win32':
progEnv.Tool('registerTargets', package = 'CRflux',
             libraryCxts = [[CRflux, libEnv]],
             testAppCxts = [[test_CRflux, progEnv]],
             includes = listFiles(['src/*.h', 'src/*.hh']),
             xml = ['xml/source_library.xml', 'xml/source_library_OpsSim.xml'],
             jo=['src/test/jobOptions.txt'])
	





