#!/bin/bash
RELEASE_VERSION=$1
wget -O src/gopeaks https://github.com/maxsonBraunLab/gopeaks/releases/download/${RELEASE_VERSION}/gopeaks-linux-amd64
chmod +x src/gopeaks
