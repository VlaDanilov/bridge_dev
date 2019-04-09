#!/bin/sh

if [ -e  src/User/gridwrap.f ];then 
  rm  src/User/gridwrap.f
else
  echo "src/User/gridwrap.f already removed"
fi

