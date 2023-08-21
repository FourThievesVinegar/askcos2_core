export VERSION ?= dev
export TAG ?= $(VERSION)

.PHONY: pre-deploy deploy start down prune

pre-deploy:
	bash deploy.sh pre-deploy -v $(TAG)

deploy:
	bash deploy.sh deploy -v $(TAG)

update:
	bash deploy.sh deploy -n -v $(TAG)

start:
	bash deploy.sh start

stop:
	bash deploy.sh stop

clean:
	bash deploy.sh clean

restart:
	bash deploy.sh restart
