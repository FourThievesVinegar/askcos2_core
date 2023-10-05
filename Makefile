export VERSION ?= dev
export TAG ?= $(VERSION)

ASKCOS_VUE_NGINX_DIR ?= ../../askcos-vue-nginx

.PHONY: pre-deploy deploy update start stop clean restart test test-wrappers test-adapters build-vue

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

test: | test-adapters test-wrappers

test-adapters:
	pytest tests/adapters

test-wrappers:
	pytest tests/wrappers

build-vue:
	$(MAKE) -C $(ASKCOS_VUE_NGINX_DIR) VERSION=$(VERSION) TAG=$(TAG)