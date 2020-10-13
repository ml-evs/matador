FROM python:3.7

WORKDIR /app

COPY . .

RUN pip install --upgrade pip && \
    pip install --no-cache-dir .[optional]
