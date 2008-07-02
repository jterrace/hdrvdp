#!/bin/bash

MAN_PAGES=`find ./src -name "*.1"`
DESTDIR='www_doc'

mkdir -p $DESTDIR

for FILE in ${MAN_PAGES}; do
    BASENAME=`basename $FILE`
    man2html -r $FILE >${DESTDIR}/man1/${BASENAME}.html
done
