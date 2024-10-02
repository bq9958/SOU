PYTHON = python3
SCRIPT = convection_diffusion.py
OUTPUT = result.log

all: run

run:
	$(PYTHON) $(SCRIPT)

clean:
	rm -f $(OUTPUT)
