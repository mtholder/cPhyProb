#!/bin/sh
set -x
d=`dirname $0`
docd="$d/doc"
cd "$docd" || exit
invocation="rst2html.py --stylesheet-path=./cPhyProb.css --link-stylesheet"
$invocation ../Authors.txt Authors.html
$invocation ../changelog.txt ChangeLog.html
$invocation ../FAQ.txt FAQ.html
$invocation ../INSTALL Installation.html
$invocation ../MAINTAINERS Maintainers.html
$invocation ../News.txt News.html
$invocation ../README index.html
$invocation ../Thanks.txt Thanks.html
$invocation ../Todo.txt Todo.html
cp ../GPL.txt ./
