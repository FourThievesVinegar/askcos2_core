export VERSION ?= dev
export TAG ?= $(VERSION)

.PHONY: pre-deploy deploy update start stop clean restart

pre-deploy:
	bash deploy.sh pre-deploy -v $(TAG)

deploy:
	bash deploy.sh deploy -v $(TAG)

update:
	bash deploy.sh update -n -v $(TAG)

start:
	bash deploy.sh start

stop:
	bash deploy.sh stop

clean:
	bash deploy.sh clean

restart:
	bash deploy.sh restart
