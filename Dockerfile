#base
FROM --platform=linux/x86_64 ubuntu:22.04

# Set the working directory inside the container
WORKDIR /opt

#run system dependencies
RUN apt-get update && apt-get install -y \
    wget unzip git g++ make python3 python3-pip \
    zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev 

# 3. Install Eagle
RUN wget https://data.broadinstitute.org/alkesgroup/Eagle/downloads/Eagle_v2.4.1.tar.gz \
    && tar -xvzf Eagle_v2.4.1.tar.gz \
    && mv Eagle_v2.4.1/eagle /usr/local/bin/ \
    && rm -rf Eagle_v2.4.1 Eagle_v2.4.1.tar.gz
    
# Install built tools 
RUN apt-get update && apt-get install -y \
    git \
    autoconf \
    automake \
    make \
    g++ \
    && rm -rf /var/lib/apt/lists/*

# 4. Install RFMix
RUN git clone https://github.com/slowkoni/rfmix.git /opt/rfmix \
    && cd /opt/rfmix \
    && autoreconf --force --install \
    && ./configure \
    && make

ENV PATH="/opt/rfmix/bin:${PATH}"

#clean
RUN apt-get clean
RUN rm -rf /var/lib/apt/lists/*

# default path
ENV PATH="/usr/local/bin:$PATH"

CMD ["bash"]
