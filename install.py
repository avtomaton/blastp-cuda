#!/usr/bin/python3

import argparse
import fileinput
import os
import shutil
import subprocess
import time


class Task(object):
    def __str__(self):
        return str(self.__dict__)

task = Task()

arg_parser = argparse.ArgumentParser(description='CUDA BLASTP installer.')
arg_parser.add_argument('--src-dir', required=True, help='sources directory')
arg_parser.add_argument('--blast-dir', default=None, help='NCBI BLAST sources directory')
arg_parser.add_argument('--blast-download', action='store_true', help='download NCBI BLAST sources')
arg_parser.add_argument('--build-dir', default='.', help='build directory')
arg_parser.add_argument('--log-dir', default='log', help='log directory')
arg_parser.add_argument('--action', choices=['install', 'diff'], default='install',
                        help='patch and build BLAST or simply view diff')

blast_options = ['--without-debug', '--with-mt', '--without-sybase', '--without-ftds',
                 '--without-fastcgi', '--without-ncbi-c', '--without-sssdb', '--without-sss',
                 '--without-geo', '--without-sp', '--without-orbacus', '--without-boost']

arg_parser.parse_args(namespace=task)

src_dir = os.path.abspath(task.src_dir)
if task.blast_dir is None:
    blast_dir = src_dir + "/../blast-ncbi-2.2.28/c++"
else:
    blast_dir = os.path.abspath(task.blast_dir)
replacement_dir = src_dir + "/ncbi_blast_files"
build_dir = os.path.abspath(task.build_dir)
log_dir = os.path.abspath(task.log_dir)

if os.path.isdir(build_dir + "/ncbi_blast_files"):
    raise RuntimeError("possibly you run this script from source dir, DO NOT DO IT!")
if not os.path.isdir(replacement_dir):
    raise RuntimeError("can't find CUDA BLASTP sources in '" + replacement_dir + "'")

if not os.path.exists(build_dir):
    os.makedirs(build_dir)
shutil.rmtree(build_dir + "/patch-src", ignore_errors=True)
shutil.rmtree(build_dir + "/ncbi-blast-src", ignore_errors=True)
shutil.rmtree(build_dir + "/bin", ignore_errors=True)
shutil.rmtree(build_dir + "/build", ignore_errors=True)
shutil.rmtree(build_dir + "/inc", ignore_errors=True)
shutil.rmtree(build_dir + "/lib", ignore_errors=True)
shutil.rmtree(build_dir + "/log", ignore_errors=True)

print('Copying sources...')
shutil.copytree(replacement_dir, build_dir + "/patch-src")

if task.blast_download:
    print('Downloading BLAST...')
    if os.path.exists(blast_dir):
        shutil.rmtree(blast_dir, ignore_errors=True)
    os.chdir(blast_dir + '/..')
    p = subprocess.Popen(['wget',
                          'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.28/ncbi-blast-2.2.28+-src.tar.gz'])
    p.wait()
    p = subprocess.Popen(['tar', '-zxf', 'ncbi-blast-2.2.28+-src.tar.gz'])
    p.wait()
    shutil.move('ncbi-blast-2.2.28+-src', build_dir + "/ncbi-blast-src")
else:
    if not os.path.isdir(blast_dir):
        raise RuntimeError("can't find vanilla BLAST sources in '" + blast_dir + "'")
    shutil.copytree(blast_dir, build_dir + "/ncbi-blast-src")

replacement_dir = build_dir + "/patch-src"
blast_dir = build_dir + "/ncbi-blast-src"
blast_options.append('--with-build-root=' + build_dir)

files_add = {'gpu_cpu_common.h': 'include/algo/blast/core', 'functions.cpp': 'src/algo/blast/core'}
files_replace = {
    'blast_engine.c': 'src/algo/blast/core',
    'blast_engine.h': 'include/algo/blast/core',
    'blast_seqsrc.c': 'src/algo/blast/core',
    'seqsrc_seqdb.cpp': 'src/algo/blast/api',
    'seqdb.cpp': 'src/objtools/blast/seqdb_reader',
    'seqdb.hpp': 'include/objtools/blast/seqdb_reader',
    'seqdbimpl.cpp': 'src/objtools/blast/seqdb_reader',
    'seqdbimpl.hpp': 'src/objtools/blast/seqdb_reader',
    'cmdline_flags.cpp': 'src/algo/blast/blastinput',
    'cmdline_flags.hpp': 'include/algo/blast/blastinput',
    'blast_args.cpp': 'src/algo/blast/blastinput',
    'blast_args.hpp': 'include/algo/blast/blastinput',
    'blastp_args.cpp': 'src/algo/blast/blastinput',
    'blast_options_cxx.cpp': 'src/algo/blast/api',
    'blast_options.hpp': 'include/algo/blast/api',
    'blast_options_local_priv.hpp': 'src/algo/blast/api',
     # 'bl2seq.cpp': 'src/algo/blast/api',
    'blast_memento_priv.hpp': 'src/algo/blast/api',
    'prelim_search_runner.hpp': 'src/algo/blast/api',
    'blast_options.h': 'include/algo/blast/core',
    'blast_seqsrc_impl.h': 'include/algo/blast/core',

    # files related to makeblastdb
    'makeblastdb.cpp': 'src/app/blastdb',
    'writedb.cpp': 'src/objtools/blast/seqdb_writer',
    'writedb.hpp': 'include/objtools/blast/seqdb_writer',
    'build_db.cpp': 'src/objtools/blast/seqdb_writer',
    'build_db.hpp': 'include/objtools/blast/seqdb_writer',
    'writedb_impl.cpp': 'src/objtools/blast/seqdb_writer',
    'writedb_impl.hpp': 'src/objtools/blast/seqdb_writer',
    'writedb_volume.cpp': 'src/objtools/blast/seqdb_writer',
    'writedb_volume.hpp': 'src/objtools/blast/seqdb_writer',
    'writedb_files.hpp': 'include/objtools/blast/seqdb_writer'
    }

for file, path in files_add.items():
    shutil.copy(replacement_dir + '/' + file, blast_dir + '/' + path + '/' + file)

print('Patching sources...')
for file, path in files_replace.items():
    new_file = replacement_dir + '/' + file
    old_file = blast_dir + '/' + path + '/' + file
    if not os.path.exists(old_file):
        raise RuntimeError(old_file + " cannot be found. Exiting")
    if task.action == 'install':
        shutil.copy(new_file, old_file)
    elif task.action == 'diff':
        p = subprocess.Popen(['kdiff3', new_file, old_file])
        p.wait()

if task.action != 'install':
    exit(0)

if not os.path.exists(log_dir):
    os.makedirs(log_dir)

print('Configuring BLAST...')
os.chdir(blast_dir)
out = open(log_dir + '/configure-plain.out', 'wb')
p = subprocess.Popen(['./configure'] + blast_options, stderr=out, stdout=out)
p.wait()
out.close()
if p.returncode != 0:
    raise RuntimeError("Error while configuring BLAST. See 'configure-plain.out'")

print('Compiling CUDA code...')
os.chdir(src_dir)
out = open(log_dir + '/gpu.out', 'wb')
p = subprocess.Popen(['make', '-f', 'Makefile', '-B',
                      'BUILD_DIR=' + build_dir, 'BLAST_DIR=' + blast_dir,
                      'PATCH_DIR=' + replacement_dir],
                     stderr=out, stdout=out)
p.wait()
out.close()
if p.returncode != 0:
    raise RuntimeError("Error while building CUDA code. See 'gpu.out'")

print('Reconfiguring BLAST for building with CUDA code...')
with fileinput.input(build_dir + '/build/Makefile.mk', inplace=True) as f:
    for line in f:
        if line.startswith('CONF_LIBS') and '-lgpublast' not in line:
            print(line.strip() + ' -lgpublast -lcuda -lcudart')
        else:
            print(line, end='')

print('Building BLAST (can be VERY long, be patient)...')
os.chdir(build_dir + '/build')
out = open(log_dir + '/build.out', 'wb')
p = subprocess.Popen(['make', 'all_r'], stderr=out, stdout=out)
while True:
    p.poll()
    print('.', end='')
    if p.returncode is not None:
        break
    time.sleep(1)
out.close()
print()
if p.returncode != 0:
    raise RuntimeError("Error while building patched BLAST code. See 'build.out'")

print("Done")
