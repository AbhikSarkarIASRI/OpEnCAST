FROM rocker/shiny:4.4.3

RUN apt-get update && apt-get install -y python3 python3-pip
RUN pip3 install streamlit

RUN R -e "install.packages(c('seqinr','ranger','gbm','e1071','xgboost','ftrCOOL'))"

COPY shiny-app /srv/shiny-server/
COPY app.py /app.py

EXPOSE 8501
CMD streamlit run /app.py
