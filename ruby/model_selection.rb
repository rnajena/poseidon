#!/usr/bin/env ruby

class ModelSelection

  attr_reader :model

  MODEL_SELECTION_BATCH = '../tools/hyphy/lib/hyphy/TemplateBatchFiles/ModelTest.bf'
  OPENMPI = '../tools/openmpi/bin/mpirun'
  HYPHY = '../tools/hyphy/bin/HYPHYMPI'

  def initialize(aln, tree, rate_classes, model_selection_method, model_rejection_level, model_out_dir, threads, parameter_string)
    @model = '010010' # the default HKY85 model, use this if something happens

    model_log_file = run_model_selection(aln, tree, rate_classes, model_selection_method, model_rejection_level, model_out_dir, threads, parameter_string)

    get_model_string(model_log_file) if model_log_file
  end

  def run_model_selection(aln, tree, rate_classes, model_selection_method, model_rejection_level, model_out_dir, threads, parameter_string)
    parameter_string << "(echo \"#{aln}\"; echo \"#{tree}\"; echo \"#{rate_classes}\"; echo \"#{model_selection_method}\"; echo \"#{model_rejection_level}\"; echo \"#{model_out_dir}/#{File.basename(aln)}.result\") | mpirun -np #{threads} HYPHYMPI ModelTest.bf\n\n"
    unless File.exists?("#{model_out_dir}/model.log")
      model_log = File.open("#{model_out_dir}/model.log",'w')
      model_log << `(echo "#{aln}"; echo "#{tree}"; echo "#{rate_classes}"; echo "#{model_selection_method}"; echo "#{model_rejection_level}"; echo "#{model_out_dir}/#{File.basename(aln)}.result") | #{OPENMPI} -np #{threads} #{HYPHY} #{MODEL_SELECTION_BATCH}`
      model_log.close
      model_log.path
    end
  end

  def get_model_string(model_log_file)
    File.open(model_log_file,'r').each do |l|
      if l.start_with?('Model String:')
        @model = l.split(':')[1].chomp
      end
    end
  end

end