#! /bin/python3

from os import system

for i in range(30):
	system("poetry run pytest")
