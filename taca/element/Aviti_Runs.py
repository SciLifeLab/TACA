from taca.element.Element_Runs import Run


class Aviti_Run(Run):
    def __init__(self, run_dir, configuration):
        self.sequencer_type = "Aviti"
        self.demux_dir = "Demultiplexing"
        super().__init__(run_dir, configuration)
