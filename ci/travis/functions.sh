#!/bin/sh

log() {
	tput bold
	tput setaf 2
	printf "%12s " $1
	tput sgr0
	shift 1
	echo $@
}

error() {
	tput bold
	tput setaf 1
	printf "%12s " $1
	tput sgr0
	shift 1
	echo $@
}
