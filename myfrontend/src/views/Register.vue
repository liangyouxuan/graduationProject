<template>
  <div class="register-container">
    <el-card class="register-card">
      <template #header>
        <h2>注册</h2>
      </template>
      
      <el-form :model="registerForm" :rules="rules" ref="registerFormRef">
        <el-form-item prop="username">
          <el-input 
            v-model="registerForm.username"
            placeholder="用户名"
            prefix-icon="User"
          />
        </el-form-item>
        
        <el-form-item prop="password">
          <el-input
            v-model="registerForm.password" 
            type="password"
            placeholder="密码"
            prefix-icon="Lock"
            show-password
          />
        </el-form-item>

        <el-form-item prop="confirmPassword">
          <el-input
            v-model="registerForm.confirmPassword"
            type="password" 
            placeholder="确认密码"
            prefix-icon="Lock"
            show-password
          />
        </el-form-item>

        <el-form-item>
          <el-button 
            @click="handleRegister" 
            class="custom-button"
            type="success"
          >
            注册
          </el-button>
        </el-form-item>

        <div class="login-link">
          已有账号? 
          <router-link to="/login">立即登录</router-link>
        </div>
      </el-form>
    </el-card>
  </div>
</template>

<script setup>
import { ref, reactive } from 'vue'
import { useRouter } from 'vue-router'
import { ElMessage } from 'element-plus'
import axios from 'axios'

const router = useRouter()
const registerFormRef = ref()

const registerForm = reactive({
  username: '',
  password: '',
  confirmPassword: ''
})

const validatePass2 = (rule, value, callback) => {
  if (value !== registerForm.password) {
    callback(new Error('两次输入密码不一致!'))
  } else {
    callback()
  }
}

const rules = {
  username: [
    { required: true, message: '请输入用户名', trigger: 'blur' },
    { min: 3, message: '用户名长度至少为3个字符', trigger: 'blur' }
  ],
  password: [
    { required: true, message: '请输入密码', trigger: 'blur' },
    { min: 6, message: '密码长度至少为6个字符', trigger: 'blur' }
  ],
  confirmPassword: [
    { required: true, message: '请再次输入密码', trigger: 'blur' },
    { validator: validatePass2, trigger: 'blur' }
  ]
}

const handleRegister = async () => {
  if (!registerFormRef.value) return

  await registerFormRef.value.validate(async (valid) => {
    if (valid) {
      try {
        const response = await axios.post('/api/register/', {
          username: registerForm.username,
          password: registerForm.password
        })
        
        if (response.status === 201) {
          ElMessage.success('注册成功')
          router.push('/login')
        }
      } catch (error) {
        let errorMessage = '注册失败'
        if (error.response?.data) {
          // 处理具体的错误信息
          if (error.response.data.username) {
            errorMessage = error.response.data.username[0]
          } else if (error.response.data.password) {
            errorMessage = error.response.data.password[0]
          } else if (typeof error.response.data === 'string') {
            errorMessage = error.response.data
          }
        }
        ElMessage.error(errorMessage)
      }
    }
  })
}
</script>

<style scoped>
.register-container {
  display: flex;
  justify-content: center;
  align-items: center;
  height: 100vh;
  background: linear-gradient(120deg, #84fab0 0%, #8fd3f4 100%);
}

.register-card {
  width: 400px;
  background: rgba(255, 255, 255, 0.9);
  border-radius: 15px;
  box-shadow: 0 8px 20px rgba(0, 0, 0, 0.1);
}

.login-link {
  text-align: center;
  margin-top: 15px;
  font-size: 14px;
}

.login-link a {
  color: #2c785c;
  text-decoration: none;
  font-weight: bold;
}

.login-link a:hover {
  color: #84fab0;
}

.custom-button {
  width: 100%;
  background-color: #2c785c !important;
  border-color: #2c785c !important;
}

.custom-button:hover {
  background-color: #84fab0 !important;
  border-color: #84fab0 !important;
  color: white !important;
}
</style> 