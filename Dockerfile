FROM python:3-slim
ENV PYTHONUNBUFFERED 1
ENV PIP_NO_CACHE_DIR 1
RUN mkdir /code
WORKDIR /code
COPY requirements.txt /code/
RUN pip install -r requirements.txt
COPY ./Codes/*.py /code/
COPY ["main.py", "/" ]
CMD ["/bin/sh"]
