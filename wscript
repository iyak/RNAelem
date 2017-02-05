from subprocess import Popen,PIPE
import sys

out="./build"

def exe(cmd):
    return Popen(cmd.split(), stdout=PIPE).stdout.read().decode("utf8").strip()

def options(opt):
    opt.load("compiler_cxx waf_unit_test")

def configure(cnf):
    cnf.load("compiler_cxx waf_unit_test")
    cnf.find_program("freetype-config", var="FTCNF")

    cnf.check_cfg(
            path="freetype-config",
            args="--cflags --libs",
            package="",
            uselib_store="freetype"
            )

    if (sys.platform.startswith("win") or sys.platform.startswith("cygwin")):
        print("for windows, I don't get path for fonts.")
        pass
    elif (sys.platform.startswith("darwin")):
        ttfs = exe("find /Library/Fonts -name *.ttf")
        for ttf in ttfs.split("\n"):
            if ("gothic" in ttf or "Gothic" in ttf):
                cnf.env.append_unique("DEFINES", ["_DFONT=\"%s\""%ttf])
                print("set font:", ttf)
                break
        else:
            print("could not find font file.")
            print("please set manually")
    else:
        ttfs = exe("find /usr/share/fonts -name *.ttf")
        for ttf in ttfs.split("\n"):
            if ("gothic" in ttf or "Gothic" in ttf):
                cnf.env.append_unique("DEFINES", ["_DFONT=\"%s\""%ttf])
                print("set font:", ttf)
                break
        else:
            print("could not find font file.")
            print("please set manually")


def build(bld):
    bld.program(
          features="cxx cxxprogram",
          cxxflags="-std=c++14 -Wall -O3",
          source="RNAelem/main.cpp",
          includes="RNAelem",
          target="bin/RNAelem",
          use="freetype")

    bld(
          features="cxx",
          cxxflags="-pthread -c",
          source="RNAelem-test/gtest/src/gtest-all.cc",
          includes="RNAelem-test/gtest"
          " RNAelem-test/gtest/include",
          target="gtest")

    bld.program(
          features="test",
          cxxflags="-std=c++14 -Wall -O3",
          linkflags="-pthread", # necessary for icpc
          source="RNAelem-test/test.cpp",
          includes="RNAelem RNAelem-test"
          " RNAelem-test/gtest/include",
          target="bin/RNAelem-test",
          use="gtest freetype")

    bld.program(
          features="test",
          cxxflags="-std=c++14 -Wall -O3",
          linkflags="-pthread", # necessary for icpc
          source="RNAelem-test/test-exact.cpp",
          includes="RNAelem RNAelem-test"
          " RNAelem-test/gtest/include",
          target="bin/RNAelem-test-exact",
          use="gtest freetype")

def test(ctx):
    ctx.exec_command("build/bin/RNAelem-test --gtest_color=yes")
    ctx.exec_command("build/bin/RNAelem-test-exact --gtest_color=yes")
