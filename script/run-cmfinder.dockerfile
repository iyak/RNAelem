FROM ubuntu:16.04

RUN ldconfig -p > /lib.before.list

RUN apt-get update && apt-get install -y \
	g++ \
	make \
    wget \
	&& rm -rf /var/lib/apt/lists/*

# can be downloaded from https://sourceforge.net/projects/weinberg-cmfinder/
WORKDIR /
# RUN wget https://sourceforge.net/projects/weinberg-cmfinder/files/cmfinder-0.4.1.18.tgz
# RUN tar -xzf cmfinder-0.4.1.18.tgz
# WORKDIR cmfinder-0.4.1.18
# RUN ./configure --enable-threads
# RUN make -j

# RUN wget http://eddylab.org/infernal/infernal-1.1.5-linux-intel-gcc.tar.gz
# RUN tar -xzf infernal-1.1.5-linux-intel-gcc.tar.gz
# RUN ldconfig -p > /lib.after.list
# RUN sort /lib.before.list > /lib.before.list.sorted
# RUN sort /lib.after.list > /lib.after.list.sorted
# RUN comm -23 /lib.after.list.sorted /lib.before.list.sorted

# RUN find . -type f -executable -exec ldd {} \; 2>/dev/null | grep "=>" | sed 's/(.*)//g' | sort | uniq

WORKDIR /
COPY CMfinder_0.2.tgz .
RUN tar -xzf CMfinder_0.2.tgz
WORKDIR /CMfinder_0.2
ENV PATH $PATH:/CMfinder_0.2/bin
ENV CMfinder /CMfinder_0.2/bin


