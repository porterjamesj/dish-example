#!/bin/bash

# add multiverse
sed -i "/^# deb.*multiverse/ s/^# //" /etc/apt/sources.list
apt-get update

# mount glusterfs
apt-get install --yes cifs-utils
mkdir /glusterfs

export TUKEY_METADATA=http://tukey-meta-data:6666/modules/v0/

SAMBA_USER=$(wget -O - -q ${TUKEY_METADATA}username)
SAMBA_PASS=$(wget -O - -q ${TUKEY_METADATA}password)

mount -t cifs \\\\cloud-controller\\glusterfs /glusterfs/ -o user=$SAMBA_USER,password=$SAMBA_PASS,noperm,nobrl,rsize=65536

# install tuxedo tools
apt-get install --yes bowtie2 tophat cufflinks

# get trimmomatic
apt-get install --yes unzip default-jre
mkdir /usr/local/java
cd /usr/local/java
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip
unzip Trimmomatic-0.32.zip


# install python tools
cd /tmp
apt-get install --yes git g++ python2.7-dev python-setuptools python-pip zlib1g-dev python-numpy
pip install pysam
git clone https://github.com/porterjamesj/dish.git
cd dish
python setup.py install
