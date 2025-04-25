clear;clc;

wrk_dir = '/home/decision_lab/work/github/frg/';
cd(wrk_dir)
stai = readtable('frg_stai.csv');

bhv_dir = [char(wrk_dir) 'bhv/'];
cd(bhv_dir)

% bhv_read error:
% 45956 bhv_read error
% 47801_postlong absent
% 37532_post bhv_read error
% 43543_post bhv_read error
% 45528_pre bhv_read error
% 47678_pre bhv_read error
% 47678_post bhv_read error
% 47801 bhv_read error

% rt issues:
% 7873_preshort, 26383_preshort, 41302, 41645_prelong, 42125_postshort,
% 42949, 43173_prelong, 44070, 45436_prelong, 47131, 47744, 48238_postshort

bhv_badids = [37532 43543 45528 45956 47801 47678]; 

bhv_fnames = dir('./*.bhv');
fname_split = split({bhv_fnames(:).name},'_');
subids = unique(str2double(fname_split(:,:,1)));
subids = setdiff(subids,bhv_badids);

stress_conditions = {'pre','post'};
tt_envs = {'short','long'};

cd(wrk_dir)

%%
% stress_conditions = {'post'};
% subids = [47801];

bhv_tab = table();
for sub_idx = 1:length(subids)
    sub_tab = table();
    subid = subids(sub_idx);
    sub_stai = table2array(stai(stai.ParticipantID==subid,'STAI_Trait'));
    for sc_idx = 1:length(stress_conditions)
        sc = stress_conditions{sc_idx};
        bhv_fname = [bhv_dir,num2str(subid),'_',sc,'stress.bhv'];
        try
            sub_bhv = bhv_read(bhv_fname);

            trial_rewards = {sub_bhv.UserVars.rwrd};
            trial_rewards(cellfun(@isempty, trial_rewards)) = {0.0};
            trial_rewards = cell2mat(trial_rewards);
            blk_leave_idxs = find(diff(trial_rewards)>5); %overstay even after rew=0
            num_patches = length(blk_leave_idxs);last_leave = max(blk_leave_idxs);

            [blk_trials,blk_patch_lens,blk_patch_trials,choices] = leave2choice(blk_leave_idxs);
            choice = cell2mat(choices);

            trial_rt = sub_bhv.ReactionTime;
            env_trial_split = min(find(isnan(trial_rt)));

            trial_times = diff(sub_bhv.AbsoluteTrialStartTime');
            travel_times = trial_times(blk_leave_idxs);
            tt_long = find(travel_times>2e4);
            tt_short = find(~(travel_times>2e4));
            if max(tt_long)<max(tt_short)
                ord = [2 1];
                
                num_patches_long = length(blk_leave_idxs(blk_leave_idxs<=env_trial_split));
                num_patches_short = length(blk_leave_idxs)-num_patches_long;
                trials_long = cell2mat(blk_trials(1:num_patches_long));
                trials_short = cell2mat(blk_trials(num_patches_long+1:end));
                num_trials_long = length(trials_long);num_trials_short = length(trials_short);

                patch_trials_long = cell2mat(blk_patch_trials(1:num_patches_long));
                patch_trials_short = cell2mat(blk_patch_trials(num_patches_long+1:end));
                patch_len_long = blk_patch_lens(1:num_patches_long);
                patch_len_short = blk_patch_lens(num_patches_long+1:end);

                times_long = trial_times(1:num_trials_long);
                times_short = trial_times(num_trials_long+1:last_leave);

                rew_long = trial_rewards(1:num_trials_long);
                rew_short = trial_rewards(num_trials_long+1:last_leave);
                rt_long = trial_rt(1:num_trials_long);
                rt_short = trial_rt(num_trials_long+1:last_leave);
                choice_long = choice(1:num_trials_long);
                choice_short = choice(num_trials_long+1:last_leave);
            else
                ord = [1 2];

                num_patches_short = length(blk_leave_idxs(blk_leave_idxs<=env_trial_split));
                num_patches_long = length(blk_leave_idxs)-num_patches_short;
                trials_short = cell2mat(blk_trials(1:num_patches_short));
                trials_long = cell2mat(blk_trials(num_patches_short+1:end));
                num_trials_long = length(trials_long);num_trials_short = length(trials_short);
                
                patch_trials_short = cell2mat(blk_patch_trials(1:num_patches_short));
                patch_trials_long = cell2mat(blk_patch_trials(num_patches_short+1:end));
                patch_len_short = blk_patch_lens(1:num_patches_short);
                patch_len_long = blk_patch_lens(num_patches_short+1:end);

                times_short = trial_times(1:num_trials_short);
                times_long = trial_times(num_trials_short+1:last_leave);

                rew_short = trial_rewards(1:num_trials_short);
                rew_long = trial_rewards(num_trials_short+1:last_leave);
                rt_short = trial_rt(1:num_trials_short);
                rt_long = trial_rt(num_trials_short+1:last_leave);
                choice_short = choice(1:num_trials_short);
                choice_long = choice(num_trials_short+1:last_leave);
            end

            nt_ord = [num_trials_short;num_trials_long];
            pl_ord = {patch_len_short;patch_len_long}; %%
            trial_ord = {trials_short;trials_long};

            patch_trial_ord = {patch_trials_short;patch_trials_long};
            ch_ord = {choice_short;choice_long};
            times_ord = {times_short;times_long};
            rt_ord = {rt_short;rt_long};
            rew_ord = {rew_short;rew_long};
        catch
            sprintf('%s',[num2str(subid),'_',sc,' bhv_read error'])
        end
        for env_idx = ord
            env = tt_envs{env_idx};
            sctt = [sc,env];

            env_num_trials = nt_ord(env_idx);
            blk_trial = trial_ord{env_idx};
            env_trial = trial_ord{env_idx}-trial_ord{env_idx}(1)+1;
            patch_trial = patch_trial_ord{env_idx};
            bhv_choice = ch_ord{env_idx};
            times = times_ord{env_idx};
            rt = rt_ord{env_idx};
            rew = rew_ord{env_idx};

            sub_col = ones(env_num_trials,1)*subid;
            sc_blk_col = ones(env_num_trials,1)*sc_idx;
            tt_env_col = ones(env_num_trials,1)*env_idx;
%             stai_col = ones(env_num_trials,1)*sub_stai;

            bhv_choice_col = bhv_choice';
            blk_trial_col = blk_trial';
            env_trial_col = env_trial';
            patch_trial_col = patch_trial';
            times_col = times';
            rt_col = rt';
            rew_col = rew';

            sub_cond_tab = table(sub_col,sc_blk_col,tt_env_col,blk_trial_col,...
                         env_trial_col,patch_trial_col,times_col,...
                         bhv_choice_col,rt_col,rew_col);
            sub_tab = [sub_tab;sub_cond_tab];
        end
    end
    bhv_tab = [bhv_tab;sub_tab];
end

function [env_trials,patch_lens,patch_trials,choices] = leave2choice(leave_indexes)
    entry_idxs = [1 leave_indexes(1:end-1)+1];
    env_trials = arrayfun(@(f,g) (f:g),entry_idxs,leave_indexes,'UniformOutput',false);
    patch_lens = cell2mat(cellfun(@length,env_trials,'UniformOutput',false));
    patch_trials = arrayfun(@(x) [1:x],patch_lens,'UniformOutput',false);
    choices = arrayfun(@(x) [ones(1,x-1) 0],patch_lens,'UniformOutput',false);
end