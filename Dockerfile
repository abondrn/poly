FROM continuumio/miniconda3

RUN conda install -c bioconda --yes minimap2 samtools diamond blast infernal

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

WORKDIR $HOME/poly

COPY go.mod go.sum .

RUN go mod download

COPY . .

RUN go test ./...

ENTRYPOINT ["go"]
CMD ["annotate/annotate.go"]