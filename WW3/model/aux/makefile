############################################################
# FC, CC & BIN must be defined when called
############################################################

BIN_LIST = $(BIN)/w3adc $(BIN)/w3list $(BIN)/w3prnt $(BIN)/w3split

all: $(BIN_LIST)

$(BIN)/w3adc: w3adc.f
	$(FC) -o $(BIN)/w3adc w3adc.f

$(BIN)/w3list: w3list.f
	$(FC) -o $(BIN)/w3list w3list.f

$(BIN)/w3prnt: w3prnt.f
	$(FC) -o $(BIN)/w3prnt w3prnt.f

$(BIN)/w3split: w3split.f
	$(FC) -o $(BIN)/w3split w3split.f

$(BIN)/w3xxx: w3xxx.c
	$(CC) -o $(BIN)/w3xxx w3xxx.c

clean:
	rm -f $(BIN_LIST)
