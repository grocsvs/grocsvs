FROM python:2

RUN apt-get update && apt-get install -y graphviz

RUN wget https://github.com/grocsvs/idba/archive/1.1.3g1.tar.gz \
    	 && tar -xf 1.1.3g1.tar.gz \
	 && cd idba-1.1.3g1 \
	 && ./build.sh \
	 && ./configure \
	 && make \
	 && mv bin/idba_ud /bin

RUN wget https://github.com/lh3/bwa/releases/download/v0.7.15/bwa-0.7.15.tar.bz2 \
    	 && tar -xf bwa-0.7.15.tar.bz2 \
    	 && cd bwa-0.7.15 \
	 && make \
	 && mv bwa /bin

RUN wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 \
    	 && tar -xf samtools-1.3.1.tar.bz2 \
	 && cd samtools-1.3.1 \
	 && make install

RUN git clone https://github.com/grocsvs/grocsvs.git \
    	&& cd grocsvs \
	&& pip install .
	
CMD grocsvs