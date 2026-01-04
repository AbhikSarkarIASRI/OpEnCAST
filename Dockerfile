FROM rocker/shiny:4.4.3

RUN apt-get update && apt-get install -y python3 python3-pip nginx
RUN pip3 install streamlit

RUN R -e "install.packages(c('seqinr','ranger','gbm','e1071','xgboost','ftrCOOL'))"

COPY ShinyAPP /srv/shiny-server/
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf
COPY app.py /app.py

EXPOSE 8501 3838

CMD shiny-server & streamlit run /app.py --server.port 8501 --server.address 0.0.0.0

