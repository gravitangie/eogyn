# Makefile for our code

# compiler
CC     = gcc
CFLAGS = -g # debug
# CFLAGS += -O3

# directories
S = $(PWD)/src
O = $(PWD)/obj

# files
SRC = $(wildcard $S/*.c)
OBJ = $(patsubst $S/%.c, $O/%.o, $(SRC))

# executable
NAME = eogyn 

# rules & prerequisites 

.PHONY : all clean

all : $(NAME)
	@echo "Code successfully compiled!"

$(NAME) : $(OBJ)
	@$(CC) $(CFLAGS) $^ -o $@ 

$O/%.o : $S/%.c
	@mkdir -p $(dir $@)
	@$(CC) $(CFLAGS) -c $< -o $@

clean :
	@rm -r $(NAME) $O
	@echo "Everything has been cleaned."