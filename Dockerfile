FROM continuumio/miniconda3

RUN conda install -c bioconda --yes minimap2 samtools diamond>=2.0.13 blast>=2.10.1 infernal

# Change your version here
ENV GO_VERSION=1.18

# For ref, see: https://github.com/gitpod-io/workspace-images/blob/61df77aad71689504112e1087bb7e26d45a43d10/chunks/lang-go/Dockerfile#L10
ENV GOPATH=$HOME/go-packages
ENV GOROOT=$HOME/go
ENV PATH=$GOROOT/bin:$GOPATH/bin:$PATH
RUN curl -fsSL https://dl.google.com/go/go${GO_VERSION}.linux-amd64.tar.gz | tar xzs

RUN apt-get update && apt-get install -y gcc && rm -rf /var/lib/apt/lists/*

WORKDIR $HOME

RUN git clone https://github.com/mmcguffi/pLannotate.git

ENV PLANNOTATE_VERSION=v1.2.0
ENV DB_LOC=https://github.com/mmcguffi/pLannotate/releases/download/${PLANNOTATE_VERSION}/BLAST_dbs.tar.gz
ENV ROOT_DIR=pLannotate/plannotate
RUN curl -L -o ${ROOT_DIR}/data/BLAST_dbs.tar.gz ${DB_LOC}
RUN tar -xzf ${ROOT_DIR}/data/BLAST_dbs.tar.gz -C ${ROOT_DIR}/data/
RUN rm ${ROOT_DIR}/data/BLAST_dbs.tar.gz

WORKDIR $HOME/poly

COPY go.mod go.sum .

RUN go mod download

COPY . .

RUN go test ./...

ENTRYPOINT ["go"]
#CMD ["annotate/annotate.go"]