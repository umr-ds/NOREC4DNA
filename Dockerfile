FROM python:3.13-slim AS builder
LABEL org.opencontainers.image.authors="Peter Michael Schwarz <peter.schwarz@uni-marburg.de>"
# uwsgi-plugin-python3
COPY . /norec4dna
WORKDIR /norec4dna

RUN apt-get update -y \
 && apt-get install --no-install-recommends -y build-essential gcc ffmpeg \
 && pip3 install -r requirements.txt --no-cache-dir \
 && python -m build . --installer pip && pip install . \
 && apt-get purge -y --auto-remove build-essential \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


# squash / reduce size
FROM scratch
COPY --from=builder / /
WORKDIR /norec4dna

#ENTRYPOINT ["python", "demo_raptor_encode.py", "Dorn", "--error_correction=reedsolomon", "--repairsymbols=6", "--asdna"]
ENTRYPOINT ["python", "unified_dna_raptor.py"]
#ENTRYPOINT ["python", "optimize_dist_gd.py"]
# OUTPUT TO parallel_RU10/