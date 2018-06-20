from subprocess import Popen,PIPE
import sys

out="./build"

def exe(cmd):
    return Popen(cmd.split(), stdout=PIPE).stdout.read().decode("utf8").strip()

def options(opt):
    opt.load("compiler_cxx compiler_c waf_unit_test python")

def configure(cnf):
    cnf.load("compiler_cxx compiler_c waf_unit_test python")
    cnf.find_program("freetype-config",var="FTCNF",mandatory=False)
    cnf.find_program("convert",mandatory=False)
    cnf.find_program("rsvg-convert",mandatory=False)
    cnf.check_python_version((3,))

    cnf.check_cfg(
            path="freetype-config",
            args="--cflags --libs",
            package="",
            uselib_store="freetype"
            )

def build(bld):
    bld(
            features="c",
            cflags="-c -O3",
            source="RNAelem/ushuffle/ushuffle.c",
            target="ushuffle")

    bld.program(
            features="cxx cxxprogram",
            cxxflags="-std=c++14 -Wall -O3 -static"
            " -Wno-unknown-pragmas -ffast-math",
            source="RNAelem/main.cpp",
            includes="RNAelem RNAelem/ushuffle",
            target="bin/RNAelem",
            use="ushuffle",
            lib="pthread")

    bld.program(
            features="cxx cxxprogram",
            cxxflags="-std=c++14 -Wall -O3"
            " -Wno-unknown-pragmas -ffast-math",
            source="RNAelem-plot/main.cpp",
            includes="RNAelem RNAelem-plot",
            target="bin/RNAelem-plot",
            lib="pthread")

    bld.program(
            features="cxx cxxprogram",
            cxxflags="-std=c++14 -Wall -O3"
            " -Wno-unknown-pragmas -ffast-math",
            source="RNAelem-logo/main.cpp",
            includes="RNAelem RNAelem-logo",
            target="bin/RNAelem-logo",
            use="freetype",
            lib="pthread")

    bld(
            features="cxx",
            cxxflags="-pthread -c",
            source="RNAelem-test/gtest/src/gtest-all.cc",
            includes="RNAelem-test/gtest"
            " RNAelem-test/gtest/include",
            target="gtest")

    bld.program(
            features="test",
            cxxflags="-std=c++14 -Wall -O3"
            " -Wno-unknown-pragmas -ffast-math",
            source="RNAelem-test/test.cpp",
            includes="RNAelem RNAelem-test"
            " RNAelem-test/gtest/include",
            target="bin/RNAelem-test",
            use="gtest ushuffle",
            lib="pthread")

    bld.program(
            features="test",
            cxxflags="-std=c++14 -Wall -O3"
            " -Wno-unknown-pragmas -ffast-math",
            source="RNAelem-test/test-exact.cpp",
            includes="RNAelem RNAelem-test"
            " RNAelem-test/gtest/include",
            target="bin/RNAelem-test-exact",
            use="gtest ushuffle",
            lib="pthread")

    if 0!=bld.is_install:
        scripts=["elem","kmer-psp.py","draw_motif.py","dishuffle.py"]
        exe("cp {scripts} {bindir}".format(
            scripts=' '.join(["script/"+f for f in scripts]),
            bindir=bld.env.BINDIR))
        print("- install utility scripts")

def test(ctx):
    ctx.exec_command("build/bin/RNAelem-test --gtest_color=yes")
    ctx.exec_command("build/bin/RNAelem-test-exact --gtest_color=yes")
