#Parent image
FROM python:3.8-slim

#Set the working directory within the container
WORKDIR /app

#Copy the entire contents of the current directory into the container at /app
COPY . .

#Install  dependencies specified in requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

#Define environment variable
ENV PYTHONUNBUFFERED 1

#Command to run the ETL pipeline
CMD ["python", "extract.py", "&&", "python", "transform.py", "&&", "python", "load.py", "&&", "python", "main.py"]