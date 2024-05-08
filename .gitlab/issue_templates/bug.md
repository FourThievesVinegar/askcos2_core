name: Bug Report
description: File a bug report.
title: "[Bug]: "
labels: ["Bug"]
projects: ["MLPDS/ASKCOSv2/askcos2_core", "MLPDS/ASKCOSv2/askcos2_deloy"]
assignees:
  - octocat
body:
  - type: markdown
    attributes:
      value: |
        Thanks for taking the time to fill out this bug report!
  - type: input
    id: contact
    attributes:
      label: Contact Details
      description: How can we get in touch with you if we need more info?
      placeholder: ex. email@example.com
    validations:
      required: false
  - type: textarea
    id: what-happened
    attributes:
      label: What happened?
      description: Also tell us, what did you expect to happen?
      placeholder: Tell us what you see!
      value: "A bug happened!"
    validations:
      required: true
  - type: dropdown
    id: version
    attributes:
      label: Version
      description: What version of our software are you running?
      options:
        - 2024-04 (Default)
        - 2024-01 (Edge)
      default: 0
    validations:
      required: true
  - type: dropdown
    id: deployment
    attributes:
      label: How did you deployed the ASKCOS instance with bug?
      multiple: true
      options:
        - Public site(not self-deployed)
        - docker-compose(full)
        - docker-compose(backend full)
        - docker-compose(backend retro)
        - docker-compose(retro)
        - Kubernetes via Helm
  - type: textarea
    id: logs
    attributes:
      label: Relevant log output
      description: Please copy and paste any relevant log output. This will be automatically formatted into code, so no need for backticks.
      render: shell
