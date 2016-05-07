#This file is part of FV-NFLlib

#Copyright(C) 2016 CryptoExperts

#This program is free software; you can redistribute it and / or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software Foundation,
#Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110 - 1301 USA

#C++ Compiler
CXX         = g++
CXXFLAGS    = -std=c++11 -funroll-loops -Ofast

#Echo function
ECHO        = echo -e

#Update the flags
LFLAGS      = -lnfllib -lmpfr -lgmpxx -lgmp 
IFLAGS      = -I.. -I/usr/local/include

#Targets
TARGETS     = Test_binary_tree Test_ec_additions Test_poly_encrypt Test_simple

all: $(TARGETS)

Test_binary_tree: Test_binary_tree.cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) $< -o $@ $(LFLAGS) 

Test_ec_additions: Test_ec_additions.cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) $< -o $@ $(LFLAGS) 

Test_poly_encrypt: Test_poly_encrypt.cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) $< -o $@ $(LFLAGS)

Test_simple: Test_simple.cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) $< -o $@ $(LFLAGS)

tests: $(TARGETS)
	@$(ECHO) "\n-------------------------\n"
	@$(ECHO) "Test binary tree"
	@$(ECHO) "\n-------------------------\n"
	@./Test_binary_tree

	@$(ECHO) "\n-------------------------\n"
	@$(ECHO) "Test homomorphic elliptic curve additions"
	@$(ECHO) "\n-------------------------\n"
	@./Test_ec_additions

	@$(ECHO) "\n-------------------------\n"
	@$(ECHO) "Test encrypt/decrypt/multiply polynomials mod (plaintextSpace, x^n+1)"
	@$(ECHO) "\n-------------------------\n"
	@./Test_poly_encrypt

	@$(ECHO) "\n-------------------------\n"
	@$(ECHO) "Test simple"
	@$(ECHO) "\n-------------------------\n"
	@./Test_simple

clean:
	-rm -rf $(TARGETS) *.dSYM
